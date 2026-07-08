"""Phase-1 evidence packet matrix for real-Earth geomorphology repair.

The packets are an offline planning/diagnostic artifact.  They bind source
ledger ids, theory claims, derived metrics, reference fixtures, generated-world
audits, optimization targets, and residual policy before production tuning.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

from aevum.diagnostics.physiographic_reference import (
    METRIC_SCHEMA,
    REFERENCE_SOURCES,
)
from aevum.diagnostics.province_reference_graph import province_reference_graph_summary
from aevum.diagnostics.real_earth_hypsometry import real_earth_hypsometry_fixture_summary
from aevum.diagnostics.reference_source_ledger import source_ledger_summary


SCHEMA = "aevum.reference_evidence_packets.v1"
PACKET_SCHEMA = "aevum.reference_evidence_packet.v1"

REQUIRED_PACKET_FIELDS: tuple[str, ...] = (
    "packet_id",
    "phase_ids",
    "source_ids",
    "raw_data_policy",
    "extraction_method",
    "theory_claims",
    "derived_metrics",
    "reference_fixture",
    "generated_world_audit",
    "optimization_targets",
    "residual_policy",
    "asset_review",
    "residual_owner_layers",
)

REQUIRED_ASSET_REVIEW: tuple[str, ...] = (
    "elevation.png",
    "terrain_provinces.png",
    "continental_detail_provinces.png",
    "ocean_depth_provinces.png",
    "plates.png",
    "crust_age.png",
    "history.png",
    "timeline.png",
    "hexmap.png",
)


def _metric_keys(*groups: str) -> tuple[str, ...]:
    keys: list[str] = []
    for group in groups:
        for spec in METRIC_SCHEMA.get(group, ()):
            key = str(spec["key"])
            if key not in keys:
                keys.append(key)
    return tuple(keys)


REFERENCE_EVIDENCE_PACKETS: tuple[dict[str, Any], ...] = (
    {
        "packet_id": "R1_global_hypsometry_planform",
        "phase_ids": ("phase_1", "phase_6"),
        "source_ids": (
            "NOAA_ETOPO_2022",
            "GEBCO_GRIDDED_BATHYMETRY",
            "NATURAL_EARTH_10M_PHYSICAL",
        ),
        "raw_data_policy": (
            "Raw ETOPO/GEBCO/Natural Earth files stay outside git; commit only "
            "extraction code and small derived metric fixtures."
        ),
        "extraction_method": (
            "Equal-area global land/ocean bins, distance-to-coast depth profiles, "
            "component/ribbon/coastline metrics, and checksumed fixture JSON."
        ),
        "theory_claims": (
            "Earth-like worlds vary, but global planform and hypsometry must stay "
            "inside broad feature-class envelopes.",
            "Nearshore depth, shelf/slope/abyss partitioning, ridges, and trenches "
            "must be separated from final color-map review.",
        ),
        "derived_metrics": _metric_keys("planform", "hypsometry", "ocean_bathymetry"),
        "reference_fixture": (
            "P77.real_earth_hypsometry_extraction",
            "P79.province_reference_graph_extraction: planform extension pending",
        ),
        "generated_world_audit": (
            "P78.generated_hypsometry_envelope",
            "P90.current_world_morphology_gap_inventory",
            "P101.planform_residual_baseline",
        ),
        "optimization_targets": (
            "land_fraction 0.25-0.33",
            "major land component count 4-14",
            "largest component share <=0.60",
            "land ribbon fraction <=0.35",
            "largest coastline complexity <=8.0",
        ),
        "residual_policy": (
            "Planform failures block promotion until production repair or a "
            "stronger reference extraction narrows the envelope."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": ("planform", "compiler_render"),
    },
    {
        "packet_id": "R2_province_crust_sediment_basement",
        "phase_ids": ("phase_1", "phase_2", "phase_3", "phase_6"),
        "source_ids": (
            "USGS_PHYSIOGRAPHIC_DIVISIONS_US",
            "NPS_PHYSIOGRAPHIC_PROVINCES",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "CRUST1_0",
            "GLIM_GLOBAL_LITHOLOGY",
            "NOAA_TOTAL_SEDIMENT_THICKNESS",
        ),
        "raw_data_policy": (
            "Raw province, crust, lithology, and sediment vectors/rasters stay "
            "outside git; derived graph and scalar fixtures are stored."
        ),
        "extraction_method": (
            "Province-class graph extraction, crust/sediment distribution by "
            "province, basement/support class sketches, and generated-world "
            "per-continent comparison."
        ),
        "theory_claims": (
            "Stable cratons and shields are old and strong, not automatically "
            "high plateaus.",
            "Covered platforms, intracratonic basins, forelands, failed rifts, "
            "and passive-margin lowlands lower or partition interiors through "
            "accommodation and sediment history.",
        ),
        "derived_metrics": _metric_keys("province_architecture", "hypsometry"),
        "reference_fixture": (
            "P79.province_reference_graph_extraction",
            "P83.crust_sediment_province_coupling",
        ),
        "generated_world_audit": (
            "P80.generated_major_continent_province_graph",
            "P96.crust_sediment_surface_ordering",
            "P101.crust_sediment_high_flat_repair",
        ),
        "optimization_targets": (
            "high-flat interior share <=0.02 of continental land",
            "positive lowland share below 500m and 1000m per large continent",
            "large continents contain multiple province classes",
            "high platforms require plateau/orogen/plume support",
        ),
        "residual_policy": (
            "Crust/sediment amplitude failures block promotion until high flats "
            "are partitioned by process-backed provinces."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": ("crust_sediment", "province_graph"),
    },
    {
        "packet_id": "R3_boundary_wilson_deeptime",
        "phase_ids": ("phase_2", "phase_3", "phase_4", "phase_6"),
        "source_ids": (
            "GPLATES",
            "PYGPLATES_EARTHBYTE",
            "EARTHBYTE_RECONSTRUCTIONS",
            "PB2002_PLATE_BOUNDARIES",
            "GEM_GLOBAL_ACTIVE_FAULTS",
        ),
        "raw_data_policy": (
            "External plate-boundary and reconstruction data stay outside git; "
            "benchmarks store boundary-type and lifecycle fixtures only."
        ),
        "extraction_method": (
            "Boundary-type length and adjacency fixtures, scripted Wilson-cycle "
            "object lifecycles, and generated boundary-object audits."
        ),
        "theory_claims": (
            "Present-day absent features can still be required because they "
            "created current continental structure.",
            "Rifts, passive margins, sutures, transforms, ridges, and old orogens "
            "must persist through archive lineage rather than current-frame labels.",
        ),
        "derived_metrics": _metric_keys("process_parentage", "ocean_bathymetry"),
        "reference_fixture": (
            "P81.boundary_process_geometry_reference",
            "P82.wilson_cycle_lifecycle_reference",
        ),
        "generated_world_audit": (
            "P95.production_boundary_lifecycle_objects",
            "P100.integrated_reaudit_and_promotion_gate",
        ),
        "optimization_targets": (
            "ridge/transform/trench/collision diversity",
            "spherical continuity across map wrap",
            "persistent basin, margin, suture, and old-orogen lineage ids",
        ),
        "residual_policy": (
            "Boundary lifecycle residuals stay cleared only if later planform "
            "repairs preserve lineage and spherical continuity."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": ("boundary_lifecycle", "planform"),
    },
    {
        "packet_id": "R4_drainage_erosion_source_to_sink",
        "phase_ids": ("phase_1", "phase_5", "phase_6"),
        "source_ids": (
            "HYDROSHEDS_HYDROBASINS_HYDRORIVERS",
            "NOAA_TOTAL_SEDIMENT_THICKNESS",
            "GLIM_GLOBAL_LITHOLOGY",
            "NOAA_ETOPO_2022",
        ),
        "raw_data_policy": (
            "Raw hydrography, lithology, and sediment grids stay outside git; "
            "small basin/source-sink metrics and synthetic fixtures are stored."
        ),
        "extraction_method": (
            "Drainage basin graphs, divide alignment, source-to-sink budget "
            "checks, old-orogen erosion decay, and generated-world basin audits."
        ),
        "theory_claims": (
            "Surface expression should respond to drainage basins and sediment "
            "routing, not isolated cells.",
            "Old orogens can lose relief while preserving inherited province "
            "boundaries and sediment source memory.",
        ),
        "derived_metrics": _metric_keys("drainage", "process_parentage"),
        "reference_fixture": (
            "P84.source_to_sink_sediment_budget",
            "P85.drainage_divide_province_alignment",
            "P86.old_orogen_erosion_decay",
        ),
        "generated_world_audit": (
            "P97.production_drainage_source_to_sink_fields",
            "P91.integrated_real_earth_morphology_promotion_audit",
        ),
        "optimization_targets": (
            "source-to-sink sediment balance closes within tolerance",
            "drainage basins remain contiguous at region scale",
            "old-orogen relief decays without losing boundary memory",
        ),
        "residual_policy": (
            "Drainage/erosion blockers remain cleared only if later planform and "
            "interior lowering preserve basin routing and sediment closure."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": ("drainage_erosion", "crust_sediment"),
    },
    {
        "packet_id": "R5_landform_margins_mountains_plateaus",
        "phase_ids": ("phase_1", "phase_3", "phase_5", "phase_6"),
        "source_ids": (
            "GMBA_MOUNTAIN_INVENTORY",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "NOAA_ETOPO_2022",
            "GEBCO_GRIDDED_BATHYMETRY",
            "PB2002_PLATE_BOUNDARIES",
        ),
        "raw_data_policy": (
            "Raw mountain, geology, and bathymetry references stay outside git; "
            "area, relief, hierarchy, and adjacency fixtures are stored."
        ),
        "extraction_method": (
            "Mountain/plateau inventory metrics, rift-to-passive-margin sequence "
            "fixtures, shelf/slope/rise/abyss ordering, and generated object audits."
        ),
        "theory_claims": (
            "Mountains can be active ranges, old subdued belts, collision "
            "plateaus, volcanic/LIP plateaus, rift shoulders, or arcs.",
            "Passive margins should form ordered rift shoulder, escarpment, "
            "coastal lowland, shelf, slope, rise, and abyss sequences.",
        ),
        "derived_metrics": _metric_keys("province_architecture", "ocean_bathymetry"),
        "reference_fixture": (
            "P87.mountain_inventory_expression",
            "P88.rift_margin_escarpment_sequence",
            "P89.plateau_area_cap_and_decay",
        ),
        "generated_world_audit": (
            "P98.production_landform_inventory_lifecycle_fields",
            "P99.production_bathymetry_margin_sequence_fields",
            "P101 current residual scale cross-check",
        ),
        "optimization_targets": (
            "mountains object-backed and finite area",
            "plateaus parented, area-limited, and decaying",
            "passive-margin lowlands tied to shelf/sediment sequence",
        ),
        "residual_policy": (
            "Landform blockers can reappear at higher resolution; every later "
            "repair must rerun multi-resolution landform audits."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": ("landform_expression", "bathymetry_margin"),
    },
    {
        "packet_id": "R6_case_study_feature_catalog",
        "phase_ids": ("phase_1", "phase_3", "phase_5", "phase_7"),
        "source_ids": (
            "AEVUM_P69_GENERATED_BASELINE",
            "USGS_PHYSIOGRAPHIC_DIVISIONS_US",
            "NPS_PHYSIOGRAPHIC_PROVINCES",
            "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
            "EARTHBYTE_RECONSTRUCTIONS",
            "GMBA_MOUNTAIN_INVENTORY",
        ),
        "raw_data_policy": (
            "Named real-Earth case studies are stored as compact feature-class "
            "sketches and graph fixtures, not raw GIS bundles or exact map targets."
        ),
        "extraction_method": (
            "North America, South America, Africa, Eurasia, India/Australia, and "
            "ocean-basin feature-class packets with parent process and adjacency "
            "expectations."
        ),
        "theory_claims": (
            "Case studies are exemplars for feature-class coverage, not maps to "
            "copy.",
            "Every required real-Earth landform family needs either generated-world "
            "expression or a named fixture with a production path.",
        ),
        "derived_metrics": _metric_keys(
            "planform",
            "province_architecture",
            "process_parentage",
            "drainage",
            "ocean_bathymetry",
        ),
        "reference_fixture": (
            "P73.real_earth_case_study_sketches",
            "P79.province_reference_graph_extraction",
        ),
        "generated_world_audit": (
            "P91.integrated_real_earth_morphology_promotion_audit",
            "P101.current_residual_attribution",
        ),
        "optimization_targets": (
            "all required feature families covered by fixture or generated world",
            "fixture-only features remain visible blockers until production path exists",
            "wet/dry variants differ without violating feature-class logic",
        ),
        "residual_policy": (
            "Promotion remains blocked if major feature families are fixture-only "
            "without a scheduled production path."
        ),
        "asset_review": REQUIRED_ASSET_REVIEW,
        "residual_owner_layers": (
            "planform",
            "province_graph",
            "crust_sediment",
            "landform_expression",
            "bathymetry_margin",
        ),
    },
)


def _normalise_packet(packet: dict[str, Any], valid_source_ids: set[str]) -> dict[str, Any]:
    source_ids = tuple(str(source_id) for source_id in packet.get("source_ids", ()))
    missing_fields = tuple(
        field
        for field in REQUIRED_PACKET_FIELDS
        if field not in packet or packet[field] in ("", (), None)
    )
    invalid_sources = tuple(
        source_id for source_id in source_ids if source_id not in valid_source_ids
    )
    derived_metrics = tuple(str(key) for key in packet.get("derived_metrics", ()))
    theory_claims = tuple(str(item) for item in packet.get("theory_claims", ()))
    reference_fixture = tuple(str(item) for item in packet.get("reference_fixture", ()))
    generated_world_audit = tuple(
        str(item) for item in packet.get("generated_world_audit", ())
    )
    optimization_targets = tuple(
        str(item) for item in packet.get("optimization_targets", ())
    )
    asset_review = tuple(str(item) for item in packet.get("asset_review", ()))
    residual_owner_layers = tuple(
        str(item) for item in packet.get("residual_owner_layers", ())
    )
    phase_ids = tuple(str(item) for item in packet.get("phase_ids", ()))
    acceptance = {
        "required_fields_present": not missing_fields,
        "source_ids_valid": not invalid_sources and bool(source_ids),
        "raw_data_policy_explicit": bool(packet.get("raw_data_policy")),
        "theory_claims_present": len(theory_claims) >= 2,
        "derived_metrics_present": len(derived_metrics) >= 4,
        "reference_track_present": bool(reference_fixture),
        "generated_world_track_present": bool(generated_world_audit),
        "optimization_targets_present": len(optimization_targets) >= 3,
        "residual_policy_present": bool(packet.get("residual_policy")),
        "asset_review_present": set(REQUIRED_ASSET_REVIEW).issubset(asset_review),
        "residual_owner_layers_present": bool(residual_owner_layers),
        "phase1_or_later_declared": bool(phase_ids) and "phase_0" not in phase_ids,
    }
    return {
        "schema": PACKET_SCHEMA,
        "packet_id": str(packet["packet_id"]),
        "phase_ids": phase_ids,
        "source_ids": source_ids,
        "invalid_source_ids": invalid_sources,
        "raw_data_policy": str(packet["raw_data_policy"]),
        "extraction_method": str(packet["extraction_method"]),
        "theory_claims": theory_claims,
        "derived_metrics": derived_metrics,
        "reference_fixture": reference_fixture,
        "generated_world_audit": generated_world_audit,
        "optimization_targets": optimization_targets,
        "residual_policy": str(packet["residual_policy"]),
        "asset_review": asset_review,
        "residual_owner_layers": residual_owner_layers,
        "missing_required_fields": missing_fields,
        "acceptance": acceptance,
    }


def reference_evidence_packet_summary() -> dict[str, Any]:
    """Return the executable Phase-1 evidence packet matrix."""

    ledger = source_ledger_summary()
    fixture = real_earth_hypsometry_fixture_summary()
    province_graph = province_reference_graph_summary()
    valid_source_ids = {
        str(entry["source_id"]) for entry in ledger.get("entries", ())
    }
    p70_source_ids = {str(source["id"]) for source in REFERENCE_SOURCES}
    packets = tuple(
        _normalise_packet(packet, valid_source_ids)
        for packet in REFERENCE_EVIDENCE_PACKETS
    )
    packet_ids = tuple(packet["packet_id"] for packet in packets)
    duplicate_packet_ids = tuple(
        packet_id
        for packet_id, count in Counter(packet_ids).items()
        if count > 1
    )
    source_ids_used = {
        source_id for packet in packets for source_id in packet["source_ids"]
    }
    metric_keys = {key for packet in packets for key in packet["derived_metrics"]}
    metric_group_keys = {
        group: {str(spec["key"]) for spec in specs}
        for group, specs in METRIC_SCHEMA.items()
    }
    covered_metric_groups = tuple(sorted(
        group
        for group, keys in metric_group_keys.items()
        if keys & metric_keys
    ))
    residual_owner_layers = tuple(sorted({
        owner
        for packet in packets
        for owner in packet["residual_owner_layers"]
    }))
    invalid_source_refs = {
        packet["packet_id"]: packet["invalid_source_ids"]
        for packet in packets
        if packet["invalid_source_ids"]
    }
    missing_fields = {
        packet["packet_id"]: packet["missing_required_fields"]
        for packet in packets
        if packet["missing_required_fields"]
    }
    incomplete_packets = tuple(
        packet["packet_id"]
        for packet in packets
        if not all(packet["acceptance"].values())
    )
    acceptance = {
        "packet_count_sufficient": len(packets) >= 6,
        "unique_packet_ids": not duplicate_packet_ids,
        "all_required_fields_present": not missing_fields,
        "all_source_ids_valid": not invalid_source_refs,
        "p76_source_ledger_ready": ledger["status"] == "source_ledger_schema_ready",
        "p77_hypsometry_fixture_ready": fixture["status"] == "hypsometry_fixture_ready",
        "p79_province_reference_graph_ready": (
            province_graph["status"] == "province_reference_graph_ready"),
        "all_packets_complete": not incomplete_packets,
        "reference_and_generated_tracks_present": all(
            packet["acceptance"]["reference_track_present"]
            and packet["acceptance"]["generated_world_track_present"]
            for packet in packets
        ),
        "raw_data_policy_explicit": all(
            packet["acceptance"]["raw_data_policy_explicit"]
            for packet in packets
        ),
        "core_metric_groups_covered": {
            "planform",
            "hypsometry",
            "province_architecture",
            "process_parentage",
            "drainage",
            "ocean_bathymetry",
        }.issubset(set(covered_metric_groups)),
        "p101_residual_owners_covered": {"planform", "crust_sediment"}.issubset(
            set(residual_owner_layers)),
        "phase1_source_coverage_expanded": len(source_ids_used) >= 12,
        "source_ids_subset_of_p76_ledger": source_ids_used.issubset(p70_source_ids),
    }
    ready = all(acceptance.values())
    return {
        "schema": SCHEMA,
        "status": (
            "reference_evidence_packets_ready"
            if ready else "reference_evidence_packets_incomplete"
        ),
        "source_ledger_status": ledger["status"],
        "hypsometry_fixture_status": fixture["status"],
        "province_reference_graph_status": province_graph["status"],
        "packet_count": len(packets),
        "source_id_count": len(source_ids_used),
        "metric_key_count": len(metric_keys),
        "covered_metric_groups": covered_metric_groups,
        "residual_owner_layers": residual_owner_layers,
        "duplicate_packet_ids": duplicate_packet_ids,
        "invalid_source_refs": invalid_source_refs,
        "missing_required_fields": missing_fields,
        "incomplete_packets": incomplete_packets,
        "packets": packets,
        "acceptance": acceptance,
    }
