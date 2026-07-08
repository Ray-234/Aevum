"""Real-Earth case-study calibration sketches for province architecture.

The sketches are lightweight, deterministic references.  They encode feature
classes, parent processes, adjacency expectations, and broad metric envelopes
without bundling large real-Earth rasters.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

from aevum.diagnostics.physiographic_reference import REFERENCE_SOURCES


SCHEMA = "aevum.earth_case_studies.v1"

REQUIRED_CASES = (
    "north_america",
    "south_america",
    "africa",
    "eurasia",
    "australia",
)

REQUIRED_FAILURE_CATEGORIES = (
    "missing_process",
    "wrong_amplitude",
    "wrong_scale",
    "wrong_adjacency",
    "compiler_rendering_mismatch",
)

REQUIRED_ENVELOPE_KEYS = (
    "province_class_count",
    "largest_internal_province_fraction",
    "basin_or_lowland_share",
    "lowland_fraction_lt500m",
    "highland_fraction_gt2500m",
    "active_orogen_or_plateau_fraction",
    "parented_highland_fraction",
)

REQUIRED_FEATURE_CLASSES = {
    "shield",
    "platform",
    "intracratonic_basin",
    "foreland_basin",
    "active_orogen",
    "old_orogen",
    "rift_system",
    "passive_margin_lowland",
}

REQUIRED_PROCESSES = {
    "cratonization",
    "platform_subsidence",
    "intracratonic_sag",
    "collision_orogeny",
    "flexural_loading",
    "orogenic_decay",
    "continental_extension",
    "passive_margin_subsidence",
}


def _province(
    province_id: str,
    label: str,
    province_class: str,
    parent_processes: tuple[str, ...],
    role: str,
    generated_analog: str,
) -> dict[str, Any]:
    return {
        "id": province_id,
        "label": label,
        "class": province_class,
        "parent_processes": parent_processes,
        "role": role,
        "generated_analog": generated_analog,
    }


def _envelope(
    *,
    province_class_count: tuple[float, float],
    largest_internal_province_fraction: tuple[float, float],
    basin_or_lowland_share: tuple[float, float],
    lowland_fraction_lt500m: tuple[float, float],
    highland_fraction_gt2500m: tuple[float, float],
    active_orogen_or_plateau_fraction: tuple[float, float],
    parented_highland_fraction: tuple[float, float],
    passive_margin_lowland_share: tuple[float, float],
    old_orogen_share: tuple[float, float],
    rift_or_extensional_basin_share: tuple[float, float],
) -> dict[str, dict[str, float]]:
    raw = {
        "province_class_count": province_class_count,
        "largest_internal_province_fraction": largest_internal_province_fraction,
        "basin_or_lowland_share": basin_or_lowland_share,
        "lowland_fraction_lt500m": lowland_fraction_lt500m,
        "highland_fraction_gt2500m": highland_fraction_gt2500m,
        "active_orogen_or_plateau_fraction": active_orogen_or_plateau_fraction,
        "parented_highland_fraction": parented_highland_fraction,
        "passive_margin_lowland_share": passive_margin_lowland_share,
        "old_orogen_share": old_orogen_share,
        "rift_or_extensional_basin_share": rift_or_extensional_basin_share,
    }
    return {
        key: {"min": float(bounds[0]), "max": float(bounds[1])}
        for key, bounds in raw.items()
    }


CASE_STUDIES: tuple[dict[str, Any], ...] = (
    {
        "id": "north_america",
        "label": "North America",
        "source_ids": (
            "NPS_PHYSIOGRAPHIC_PROVINCES",
            "USGS_PHYSIOGRAPHIC_DIVISIONS_US",
            "NOAA_ETOPO_2022",
            "CRUST1_0",
            "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
        ),
        "provinces": (
            _province("canadian_shield", "Canadian Shield", "shield",
                      ("cratonization",), "old stable core", "shield"),
            _province("interior_plains", "Interior Plains", "platform",
                      ("platform_subsidence",), "covered stable platform", "platform"),
            _province("appalachians", "Appalachians", "old_orogen",
                      ("collision_orogeny", "orogenic_decay"), "eroded inherited orogen", "old_subdued_orogen"),
            _province("cordillera", "Cordillera", "active_orogen",
                      ("collision_orogeny", "subduction_accretion"), "active western highlands", "orogen"),
            _province("basin_and_range", "Basin and Range", "rift_system",
                      ("continental_extension",), "extensional province", "rift_basin"),
            _province("coastal_plain", "Coastal Plain", "passive_margin_lowland",
                      ("passive_margin_subsidence", "shelf_sedimentation"), "low passive margin plain", "passive_margin_lowland"),
        ),
        "adjacency": (
            ("canadian_shield", "interior_plains"),
            ("interior_plains", "appalachians"),
            ("interior_plains", "cordillera"),
            ("cordillera", "basin_and_range"),
            ("appalachians", "coastal_plain"),
        ),
        "metric_envelope": _envelope(
            province_class_count=(5, 8),
            largest_internal_province_fraction=(0.18, 0.45),
            basin_or_lowland_share=(0.25, 0.70),
            lowland_fraction_lt500m=(0.25, 0.65),
            highland_fraction_gt2500m=(0.02, 0.20),
            active_orogen_or_plateau_fraction=(0.06, 0.25),
            parented_highland_fraction=(0.85, 1.0),
            passive_margin_lowland_share=(0.04, 0.22),
            old_orogen_share=(0.04, 0.22),
            rift_or_extensional_basin_share=(0.03, 0.18),
        ),
    },
    {
        "id": "south_america",
        "label": "South America",
        "source_ids": (
            "NOAA_ETOPO_2022",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "PB2002_PLATE_BOUNDARIES",
            "CRUST1_0",
            "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
        ),
        "provinces": (
            _province("andes", "Andes", "active_orogen",
                      ("collision_orogeny", "subduction_accretion"), "active convergent-margin highlands", "orogen"),
            _province("amazon_foreland", "Amazon foreland and basin", "foreland_basin",
                      ("flexural_loading", "source_to_sink_sedimentation"), "large lowland foreland basin", "foreland_basin"),
            _province("brazilian_shield", "Brazilian Shield", "shield",
                      ("cratonization",), "old stable core", "shield"),
            _province("guiana_shield", "Guiana Shield", "shield",
                      ("cratonization",), "old northern stable core", "shield"),
            _province("patagonia", "Patagonia", "platform",
                      ("platform_subsidence", "orogenic_decay"), "southern platform and old uplands", "platform"),
            _province("atlantic_passive_margin", "Atlantic passive margin", "passive_margin_lowland",
                      ("passive_margin_subsidence", "shelf_sedimentation"), "eastern coastal lowlands and shelf", "passive_margin_lowland"),
        ),
        "adjacency": (
            ("andes", "amazon_foreland"),
            ("amazon_foreland", "brazilian_shield"),
            ("amazon_foreland", "guiana_shield"),
            ("brazilian_shield", "atlantic_passive_margin"),
            ("brazilian_shield", "patagonia"),
        ),
        "metric_envelope": _envelope(
            province_class_count=(4, 7),
            largest_internal_province_fraction=(0.18, 0.48),
            basin_or_lowland_share=(0.30, 0.72),
            lowland_fraction_lt500m=(0.30, 0.70),
            highland_fraction_gt2500m=(0.03, 0.22),
            active_orogen_or_plateau_fraction=(0.06, 0.25),
            parented_highland_fraction=(0.88, 1.0),
            passive_margin_lowland_share=(0.03, 0.20),
            old_orogen_share=(0.00, 0.16),
            rift_or_extensional_basin_share=(0.00, 0.12),
        ),
    },
    {
        "id": "africa",
        "label": "Africa",
        "source_ids": (
            "NOAA_ETOPO_2022",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "GLIM_GLOBAL_LITHOLOGY",
            "CRUST1_0",
            "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
        ),
        "provinces": (
            _province("west_central_southern_cratons", "West, Central, and Southern cratons", "shield",
                      ("cratonization",), "multiple old continental nuclei", "shield"),
            _province("congo_basin", "Congo Basin", "intracratonic_basin",
                      ("intracratonic_sag", "source_to_sink_sedimentation"), "broad low interior basin", "interior_basin"),
            _province("sahara_platform", "Sahara platform", "platform",
                      ("platform_subsidence",), "wide stable platform", "platform"),
            _province("east_african_rift", "East African Rift", "rift_system",
                      ("continental_extension", "rift_shoulder_uplift"), "active continental rift", "rift_basin"),
            _province("ethiopian_highlands", "Ethiopian Highlands", "volcanic_lip_plateau",
                      ("plume_lip_emplacement",), "volcanic plateau", "plateau"),
            _province("atlantic_indian_margins", "Atlantic and Indian passive margins", "passive_margin_lowland",
                      ("passive_margin_subsidence", "shelf_sedimentation"), "coastal and shelf lowlands", "passive_margin_lowland"),
        ),
        "adjacency": (
            ("west_central_southern_cratons", "congo_basin"),
            ("west_central_southern_cratons", "sahara_platform"),
            ("congo_basin", "east_african_rift"),
            ("east_african_rift", "ethiopian_highlands"),
            ("west_central_southern_cratons", "atlantic_indian_margins"),
        ),
        "metric_envelope": _envelope(
            province_class_count=(5, 8),
            largest_internal_province_fraction=(0.16, 0.45),
            basin_or_lowland_share=(0.25, 0.65),
            lowland_fraction_lt500m=(0.18, 0.55),
            highland_fraction_gt2500m=(0.01, 0.16),
            active_orogen_or_plateau_fraction=(0.02, 0.20),
            parented_highland_fraction=(0.85, 1.0),
            passive_margin_lowland_share=(0.04, 0.25),
            old_orogen_share=(0.02, 0.22),
            rift_or_extensional_basin_share=(0.03, 0.18),
        ),
    },
    {
        "id": "eurasia",
        "label": "Eurasia",
        "source_ids": (
            "NOAA_ETOPO_2022",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "PB2002_PLATE_BOUNDARIES",
            "GMBA_MOUNTAIN_INVENTORY",
            "CRUST1_0",
        ),
        "provinces": (
            _province("east_european_platform", "East European Platform", "platform",
                      ("platform_subsidence",), "stable covered craton", "platform"),
            _province("west_siberian_basin", "West Siberian Basin", "intracratonic_basin",
                      ("intracratonic_sag", "source_to_sink_sedimentation"), "giant low sedimentary basin", "interior_basin"),
            _province("alps_himalaya_belt", "Alps-Himalaya belt", "active_orogen",
                      ("collision_orogeny", "flexural_loading"), "active collisional mountain belt", "orogen"),
            _province("tibetan_plateau", "Tibetan Plateau", "active_orogen",
                      ("collision_orogeny", "crustal_thickening"), "large active plateau", "plateau"),
            _province("central_asian_old_orogens", "Central Asian old orogens", "old_orogen",
                      ("collision_orogeny", "orogenic_decay"), "reworked inherited belts", "old_subdued_orogen"),
            _province("arctic_and_indian_margins", "Arctic and Indian margins", "passive_margin_lowland",
                      ("passive_margin_subsidence", "shelf_sedimentation"), "marginal lowlands and shelves", "passive_margin_lowland"),
        ),
        "adjacency": (
            ("east_european_platform", "west_siberian_basin"),
            ("west_siberian_basin", "central_asian_old_orogens"),
            ("central_asian_old_orogens", "tibetan_plateau"),
            ("tibetan_plateau", "alps_himalaya_belt"),
            ("east_european_platform", "arctic_and_indian_margins"),
        ),
        "metric_envelope": _envelope(
            province_class_count=(5, 9),
            largest_internal_province_fraction=(0.14, 0.42),
            basin_or_lowland_share=(0.20, 0.62),
            lowland_fraction_lt500m=(0.20, 0.60),
            highland_fraction_gt2500m=(0.04, 0.28),
            active_orogen_or_plateau_fraction=(0.08, 0.32),
            parented_highland_fraction=(0.90, 1.0),
            passive_margin_lowland_share=(0.03, 0.18),
            old_orogen_share=(0.05, 0.25),
            rift_or_extensional_basin_share=(0.00, 0.14),
        ),
    },
    {
        "id": "australia",
        "label": "Australia",
        "source_ids": (
            "NOAA_ETOPO_2022",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "GLIM_GLOBAL_LITHOLOGY",
            "CRUST1_0",
            "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
        ),
        "provinces": (
            _province("western_shield", "Western Shield", "shield",
                      ("cratonization",), "old stable shield", "shield"),
            _province("central_lowlands", "Central Lowlands", "intracratonic_basin",
                      ("intracratonic_sag", "source_to_sink_sedimentation"), "broad low interior basin", "interior_basin"),
            _province("eastern_highlands", "Eastern Highlands", "old_orogen",
                      ("orogenic_decay", "passive_margin_uplift"), "subdued inherited highlands", "old_subdued_orogen"),
            _province("northern_platform", "Northern platform", "platform",
                      ("platform_subsidence",), "stable platform cover", "platform"),
            _province("passive_margins", "Passive margins", "passive_margin_lowland",
                      ("passive_margin_subsidence", "shelf_sedimentation"), "coastal lowlands and shelves", "passive_margin_lowland"),
        ),
        "adjacency": (
            ("western_shield", "central_lowlands"),
            ("central_lowlands", "eastern_highlands"),
            ("central_lowlands", "northern_platform"),
            ("western_shield", "passive_margins"),
            ("eastern_highlands", "passive_margins"),
        ),
        "metric_envelope": _envelope(
            province_class_count=(4, 7),
            largest_internal_province_fraction=(0.18, 0.50),
            basin_or_lowland_share=(0.35, 0.80),
            lowland_fraction_lt500m=(0.35, 0.78),
            highland_fraction_gt2500m=(0.00, 0.05),
            active_orogen_or_plateau_fraction=(0.00, 0.08),
            parented_highland_fraction=(0.80, 1.0),
            passive_margin_lowland_share=(0.05, 0.30),
            old_orogen_share=(0.03, 0.18),
            rift_or_extensional_basin_share=(0.00, 0.10),
        ),
    },
)


def earth_case_study_calibration_summary() -> dict[str, Any]:
    source_ids = {str(source["id"]) for source in REFERENCE_SOURCES}
    cases = [_analyze_case(case, source_ids) for case in CASE_STUDIES]
    all_classes = sorted({
        str(province["class"])
        for case in CASE_STUDIES
        for province in case["provinces"]
    })
    all_processes = sorted({
        str(process)
        for case in CASE_STUDIES
        for province in case["provinces"]
        for process in province["parent_processes"]
    })
    failure_categories = sorted({
        category
        for case in cases
        for category in case["failure_categories"]
    })
    acceptance = {
        "required_case_suite_complete": tuple(case["id"] for case in CASE_STUDIES) == REQUIRED_CASES,
        "all_cases_have_reference_sketch": all(case["acceptance"]["has_reference_sketch"] for case in cases),
        "all_cases_have_metric_envelope": all(case["acceptance"]["has_metric_envelope"] for case in cases),
        "all_case_sources_known": all(case["acceptance"]["source_ids_known"] for case in cases),
        "all_adjacencies_valid": all(case["acceptance"]["adjacency_edges_valid"] for case in cases),
        "feature_class_not_exact_geography_policy": all(
            case["acceptance"]["feature_class_not_exact_geography_policy"] for case in cases),
        "failure_categories_complete": all(
            case["acceptance"]["failure_categories_complete"] for case in cases),
        "required_feature_classes_covered": REQUIRED_FEATURE_CLASSES.issubset(all_classes),
        "required_parent_processes_covered": REQUIRED_PROCESSES.issubset(all_processes),
    }
    return {
        "schema": SCHEMA,
        "status": "case_study_calibration_ready" if all(acceptance.values()) else "case_study_calibration_incomplete",
        "case_count": len(cases),
        "province_count": sum(case["province_count"] for case in cases),
        "province_class_count": len(all_classes),
        "parent_process_count": len(all_processes),
        "adjacency_edge_count": sum(case["adjacency_edge_count"] for case in cases),
        "classes_covered": all_classes,
        "parent_processes_covered": all_processes,
        "failure_categories": failure_categories,
        "cases": cases,
        "acceptance": acceptance,
        "next_gates": (
            "P74.terrain_coupling_rewrite",
            "P75.release_promotion_audit",
        ),
    }


def _analyze_case(case: dict[str, Any], known_source_ids: set[str]) -> dict[str, Any]:
    province_ids = {str(province["id"]) for province in case["provinces"]}
    classes = sorted({str(province["class"]) for province in case["provinces"]})
    processes = sorted({
        str(process)
        for province in case["provinces"]
        for process in province["parent_processes"]
    })
    missing_source_ids = sorted(set(case["source_ids"]) - known_source_ids)
    invalid_edges = sorted(
        tuple(edge)
        for edge in case["adjacency"]
        if str(edge[0]) not in province_ids or str(edge[1]) not in province_ids)
    envelope = case["metric_envelope"]
    invalid_envelope_keys = [
        key for key, bounds in envelope.items()
        if float(bounds["min"]) > float(bounds["max"])
    ]
    missing_envelope_keys = sorted(set(REQUIRED_ENVELOPE_KEYS) - set(envelope))
    reproduction_policy = {
        "feature_class_required": True,
        "exact_geography_required": False,
        "metric_envelope_required": True,
    }
    failure_modes = {
        "missing_process": "Required parent process or province class is absent.",
        "wrong_amplitude": "Feature exists but elevation, relief, or sediment amplitude is outside envelope.",
        "wrong_scale": "Feature exists but occupies an implausible share of the continent.",
        "wrong_adjacency": "Feature exists but lacks the expected neighboring province relationship.",
        "compiler_rendering_mismatch": "Simulation fields pass but rendered or compiled map loses the feature.",
    }
    acceptance = {
        "has_reference_sketch": len(case["provinces"]) >= 5 and len(classes) >= 4,
        "has_metric_envelope": not missing_envelope_keys and not invalid_envelope_keys,
        "source_ids_known": not missing_source_ids,
        "adjacency_edges_valid": not invalid_edges and len(case["adjacency"]) >= 4,
        "feature_class_not_exact_geography_policy": (
            reproduction_policy["feature_class_required"]
            and not reproduction_policy["exact_geography_required"]
        ),
        "failure_categories_complete": set(REQUIRED_FAILURE_CATEGORIES).issubset(failure_modes),
        "parent_process_links_present": all(
            bool(province["parent_processes"]) for province in case["provinces"]),
    }
    return {
        "id": str(case["id"]),
        "label": str(case["label"]),
        "source_ids": tuple(str(source_id) for source_id in case["source_ids"]),
        "missing_source_ids": missing_source_ids,
        "province_count": len(case["provinces"]),
        "province_class_count": len(classes),
        "parent_process_count": len(processes),
        "adjacency_edge_count": len(case["adjacency"]),
        "classes": classes,
        "parent_processes": processes,
        "provinces": case["provinces"],
        "adjacency": case["adjacency"],
        "invalid_adjacency_edges": invalid_edges,
        "metric_envelope": envelope,
        "missing_envelope_keys": missing_envelope_keys,
        "invalid_envelope_keys": invalid_envelope_keys,
        "reproduction_policy": reproduction_policy,
        "failure_modes": failure_modes,
        "failure_categories": tuple(failure_modes),
        "acceptance": acceptance,
    }
