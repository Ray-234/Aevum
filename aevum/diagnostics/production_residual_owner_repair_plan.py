"""P92 production residual owner repair plan."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


SCHEMA = "aevum.production_residual_owner_repair_plan.v1"

PACKET_SEQUENCE = (
    "P92.1_planform_and_reference_calibration",
    "P92.2_production_province_graph_fields",
    "P92.3_boundary_lifecycle_objects",
    "P92.4_crust_sediment_interior_relief_coupling",
    "P92.5_drainage_source_to_sink_fields",
    "P92.6_landform_inventory_lifecycle",
    "P92.7_bathymetry_margin_sequence",
    "P92.8_integrated_reaudit_and_promotion_gate",
)

OWNER_TO_PACKET = {
    "planform": "P92.1_planform_and_reference_calibration",
    "province_graph": "P92.2_production_province_graph_fields",
    "boundary_lifecycle": "P92.3_boundary_lifecycle_objects",
    "crust_sediment": "P92.4_crust_sediment_interior_relief_coupling",
    "drainage_erosion": "P92.5_drainage_source_to_sink_fields",
    "landform_expression": "P92.6_landform_inventory_lifecycle",
    "bathymetry_margin": "P92.7_bathymetry_margin_sequence",
}

BLOCKER_TO_PACKET = {
    "p69_earthlike_reference_needs_calibration": (
        "P92.1_planform_and_reference_calibration",
        "P92.8_integrated_reaudit_and_promotion_gate",
    ),
    "p90_current_world_residuals_unresolved": (
        "P92.8_integrated_reaudit_and_promotion_gate",
    ),
    "planform_residuals_unresolved": (
        "P92.1_planform_and_reference_calibration",
    ),
    "province_graph_residuals_unresolved": (
        "P92.2_production_province_graph_fields",
    ),
    "boundary_lifecycle_residuals_unresolved": (
        "P92.3_boundary_lifecycle_objects",
    ),
    "crust_sediment_residuals_unresolved": (
        "P92.4_crust_sediment_interior_relief_coupling",
    ),
    "drainage_erosion_residuals_unresolved": (
        "P92.5_drainage_source_to_sink_fields",
    ),
    "landform_expression_residuals_unresolved": (
        "P92.6_landform_inventory_lifecycle",
    ),
    "bathymetry_margin_residuals_unresolved": (
        "P92.7_bathymetry_margin_sequence",
    ),
}


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def latest_p91_summary(root: Path) -> dict[str, Any] | None:
    pattern = "out_bench_p91_integrated_real_earth_morphology_promotion_audit_*/tectonics_bench_summary.json"
    candidates = sorted(root.glob(pattern))
    loaded: list[dict[str, Any]] = []
    for path in candidates:
        try:
            summary = json.loads(path.read_text())
        except Exception:
            continue
        if summary.get("schema") == "aevum.tectonics_bench.p91.v1":
            loaded.append({"path": str(path), "summary": summary})
    passing = [
        item for item in loaded
        if item["summary"].get("status") == "pass"
    ]
    return (passing or loaded)[-1] if loaded else None


def _owner_for_residual_item(item: str) -> str:
    text = str(item)
    if text == "transform" or text == "tectonics.spreading_centers":
        return "boundary_lifecycle"
    if (
        text.startswith("tectonics.continental_province")
        or text == "tectonics.province_parent_process"
    ):
        return "province_graph"
    if (
        text.startswith("terrain.drainage")
        or text.startswith("terrain.flow_")
        or text.startswith("terrain.sediment_")
        or text.startswith("terrain.old_orogen_")
        or text.startswith("terrain.orogen_")
    ):
        return "drainage_erosion"
    if (
        text.startswith("terrain.mountain")
        or text.startswith("tectonics.mountain")
        or text.startswith("terrain.plateau")
        or text.startswith("tectonics.plateau")
        or text in {"orogen", "plateau"}
    ):
        return "landform_expression"
    if (
        text.startswith("terrain.rift_margin")
        or text.startswith("tectonics.rift_margin")
        or text in {"terrain.rift_shoulders", "terrain.escarpments"}
    ):
        return "bathymetry_margin"
    return "unassigned"


def _items_for_owner(items: tuple[str, ...], owner: str) -> tuple[str, ...]:
    return tuple(item for item in items if _owner_for_residual_item(item) == owner)


def _packet_blueprints(
    *,
    owner_counts_900: dict[str, int],
    owner_counts_2500: dict[str, int],
    category_counts: dict[str, int],
    residual_items: tuple[str, ...],
    blockers: tuple[str, ...],
) -> list[dict[str, Any]]:
    residual_by_owner = {
        owner: _items_for_owner(residual_items, owner)
        for owner in OWNER_TO_PACKET
    }

    def common(owner: str) -> dict[str, Any]:
        return {
            "owner_layer": owner,
            "gap_count_900": int(owner_counts_900.get(owner, 0)),
            "gap_count_2500": int(owner_counts_2500.get(owner, 0)),
            "residual_items": residual_by_owner.get(owner, ()),
            "trigger_blockers": tuple(
                blocker for blocker in blockers
                if blocker in BLOCKER_TO_PACKET
                and OWNER_TO_PACKET[owner] in BLOCKER_TO_PACKET[blocker]
            ),
        }

    return [
        {
            "packet_id": "P92.1_planform_and_reference_calibration",
            "priority": 1,
            "status": "planned",
            **common("planform"),
            "depends_on": (),
            "implementation_targets": (
                "aevum/modules/tectonics.py planform component/ribbon/coastline controls",
                "aevum/diagnostics/earth_reference.py generated envelope comparison",
                "aevum/diagnostics/planform_reference.py planform guardrails",
            ),
            "microbenchmarks_to_add": (
                "P93.planform_reference_calibration",
                "P93.generated_component_ribbon_envelope",
            ),
            "validation_suites": ("P78", "P90", "P91"),
            "acceptance_targets": (
                "current land_fraction and component/ribbon/coastline residuals shrink without hiding p69 calibration blocker",
                "900 and 2500 cell worlds keep at least two major components without excess narrow ribbons",
                "P91 still records any remaining planform blocker by explicit metric name",
            ),
            "regression_risks": (
                "over-simplifying coastlines into artificial blobs",
                "removing valid island arcs or sutured ribbons",
            ),
        },
        {
            "packet_id": "P92.2_production_province_graph_fields",
            "priority": 2,
            "status": "planned",
            **common("province_graph"),
            "depends_on": ("P92.1_planform_and_reference_calibration",),
            "implementation_targets": (
                "aevum/modules/tectonics.py persistent province object graph",
                "terrain.continental_province_id",
                "terrain.continental_province_code",
                "tectonics.province_parent_process",
            ),
            "microbenchmarks_to_add": (
                "P94.production_province_graph_fields",
                "P94.volcanic_lip_plateau_edge_coverage",
            ),
            "validation_suites": ("P79", "P80", "P83", "P90", "P91"),
            "acceptance_targets": (
                "production worlds expose first-class province IDs and parent processes",
                "volcanic_lip_plateau and required rift/LIP edge are expressible",
                "major continents remain multi-province at 900 and 2500 cells",
            ),
            "regression_risks": (
                "duplicating raster classes without persistent province identity",
                "creating checkerboard provinces to satisfy class counts",
            ),
        },
        {
            "packet_id": "P92.3_boundary_lifecycle_objects",
            "priority": 3,
            "status": "planned",
            **common("boundary_lifecycle"),
            "depends_on": ("P92.2_production_province_graph_fields",),
            "implementation_targets": (
                "aevum/modules/tectonics.py transform process masks",
                "tectonics.spreading_centers object set",
                "archive.wilson_cycle_phase parent links",
            ),
            "microbenchmarks_to_add": (
                "P95.transform_and_spreading_center_objects",
                "P95.boundary_lifecycle_current_world_audit",
            ),
            "validation_suites": ("P81", "P82", "P90", "P91"),
            "acceptance_targets": (
                "transform boundaries appear as generated process geometry",
                "spreading-center objects persist through Wilson-cycle frames",
                "boundary changes stay spherical and seam-consistent",
            ),
            "regression_risks": (
                "turning every oblique contact into a transform",
                "inflating ridge/trench fractions while fixing object coverage",
            ),
        },
        {
            "packet_id": "P92.4_crust_sediment_interior_relief_coupling",
            "priority": 4,
            "status": "planned",
            **common("crust_sediment"),
            "depends_on": (
                "P92.2_production_province_graph_fields",
                "P92.3_boundary_lifecycle_objects",
            ),
            "implementation_targets": (
                "aevum/modules/terrain.py province-driven interior relief",
                "crust.thickness_m and sediment.thickness_m province coupling",
                "terrain.elevation_m high-flat interior reduction",
            ),
            "microbenchmarks_to_add": (
                "P96.high_flat_interior_owner_reduction",
                "P96.province_crust_sediment_surface_ordering",
            ),
            "validation_suites": ("P83", "P86", "P90", "P91"),
            "acceptance_targets": (
                "high-flat interior share drops for process reasons, not by noise",
                "basins and passive margins remain low where accommodation exists",
                "old stable crust does not become uniformly high tableland",
            ),
            "regression_risks": (
                "reintroducing fragmented line mountains",
                "erasing shield/platform distinctions while lowering interiors",
            ),
        },
        {
            "packet_id": "P92.5_drainage_source_to_sink_fields",
            "priority": 5,
            "status": "planned",
            **common("drainage_erosion"),
            "depends_on": ("P92.4_crust_sediment_interior_relief_coupling",),
            "implementation_targets": (
                "terrain.drainage_basins",
                "terrain.drainage_divides",
                "terrain.flow_direction",
                "terrain.flow_accumulation",
                "terrain.sediment_routing_edges",
                "terrain.sediment_budget",
                "terrain.old_orogen_decay_stage",
                "terrain.orogen_erosion_budget",
            ),
            "microbenchmarks_to_add": (
                "P97.production_drainage_source_to_sink_fields",
                "P97.old_orogen_decay_budget_current_world",
            ),
            "validation_suites": ("P84", "P85", "P86", "P90", "P91"),
            "acceptance_targets": (
                "drainage basins are contiguous and province-aligned",
                "flow fields reach coherent sinks without crossing divides",
                "old orogens decay by erosion/sediment budget rather than disappearing",
            ),
            "regression_risks": (
                "creating drainage grids that look coherent but ignore spherical adjacency",
                "moving too much sediment and changing land/sea masks",
            ),
        },
        {
            "packet_id": "P92.6_landform_inventory_lifecycle",
            "priority": 6,
            "status": "planned",
            **common("landform_expression"),
            "depends_on": (
                "P92.2_production_province_graph_fields",
                "P92.5_drainage_source_to_sink_fields",
            ),
            "implementation_targets": (
                "terrain.mountain_ranges",
                "terrain.mountain_inventory",
                "terrain.mountain_hierarchy_level",
                "tectonics.mountain_belt_id",
                "terrain.plateau_inventory",
                "terrain.plateau_age_myr",
                "terrain.plateau_decay_stage",
                "terrain.plateau_parent_process_id",
                "tectonics.plateau_lineage_id",
            ),
            "microbenchmarks_to_add": (
                "P98.production_mountain_inventory_fields",
                "P98.production_plateau_lifecycle_fields",
            ),
            "validation_suites": ("P87", "P89", "P90", "P91"),
            "acceptance_targets": (
                "mountain and plateau classes are backed by inventory fields",
                "plateaus remain finite and decay through time",
                "elongated mountain ranges emerge from parent processes, not thresholds alone",
            ),
            "regression_risks": (
                "overpainting all high interiors as plateaus",
                "turning every boundary into a mountain range",
            ),
        },
        {
            "packet_id": "P92.7_bathymetry_margin_sequence",
            "priority": 7,
            "status": "planned",
            **common("bathymetry_margin"),
            "depends_on": (
                "P92.3_boundary_lifecycle_objects",
                "P92.5_drainage_source_to_sink_fields",
            ),
            "implementation_targets": (
                "terrain.rift_shoulders",
                "terrain.escarpments",
                "terrain.rift_margin_sequence_id",
                "terrain.rift_margin_stage",
                "tectonics.rift_margin_lineage_id",
                "ocean.depth_province margin transitions",
            ),
            "microbenchmarks_to_add": (
                "P99.production_rift_margin_sequence_fields",
                "P99.nearshore_bathymetry_margin_profile",
            ),
            "validation_suites": ("P81", "P88", "P90", "P91"),
            "acceptance_targets": (
                "rift shoulders, escarpments, shelves, slopes, rises, and abyssal plains form ordered sequences",
                "passive-margin lowlands are not tiny artifacts",
                "nearshore superdeep anomalies remain low without making far oceans too shallow",
            ),
            "regression_risks": (
                "broadening shelves until oceans lose abyssal area",
                "misclassifying active trenches as passive margin escarpments",
            ),
        },
        {
            "packet_id": "P92.8_integrated_reaudit_and_promotion_gate",
            "priority": 8,
            "status": "planned",
            "owner_layer": "integrated_promotion",
            "gap_count_900": int(sum(owner_counts_900.values())),
            "gap_count_2500": int(sum(owner_counts_2500.values())),
            "residual_items": residual_items,
            "trigger_blockers": tuple(
                blocker for blocker in blockers
                if blocker in {
                    "p69_earthlike_reference_needs_calibration",
                    "p90_current_world_residuals_unresolved",
                }
            ),
            "depends_on": tuple(PACKET_SEQUENCE[:7]),
            "implementation_targets": (
                "aevum/diagnostics/tectonics_bench.py P91/P92 promotion gates",
                "out_bench_p91_integrated_real_earth_morphology_promotion_audit_* evidence",
            ),
            "microbenchmarks_to_add": (
                "P100.integrated_owner_repair_reaudit",
                "P100.default_promotion_decision_gate",
            ),
            "validation_suites": (
                "P78", "P80", "P81", "P82", "P83", "P84", "P85",
                "P86", "P87", "P88", "P89", "P90", "P91",
            ),
            "acceptance_targets": (
                "P90 non-asset residuals are eliminated or reduced to explicitly waived future work",
                "P91 promotion blockers drop to zero before default Earth-like promotion",
                "8000-cell P69 reference calibration no longer reports earthlike_reference_needs_calibration",
            ),
            "regression_risks": (
                "passing narrow fixture suites while current 900/2500 worlds regress",
                "promoting default Earth-like behavior before PNG and metric evidence agree",
            ),
            "category_counts": dict(category_counts),
        },
    ]


def _dependency_order_valid(packets: list[dict[str, Any]]) -> bool:
    priority = {packet["packet_id"]: int(packet["priority"]) for packet in packets}
    return all(
        dep in priority and priority[dep] < int(packet["priority"])
        for packet in packets
        for dep in packet.get("depends_on", ())
    )


def _unassigned_residual_items(items: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(item for item in items if _owner_for_residual_item(item) == "unassigned")


def production_residual_owner_repair_plan_summary(
    root: Path,
    *,
    p91_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    selected_path = ""
    if p91_summary is None:
        selected = latest_p91_summary(root)
        if selected is not None:
            selected_path = str(selected["path"])
            p91_summary = selected["summary"]
    p91_summary = p91_summary or {}
    audit = p91_summary.get("integrated_real_earth_morphology_promotion_audit", {})
    decision = audit.get("promotion_decision", {})
    ci_worlds = list(audit.get("ci_worlds", ()))
    root_owner_counts = {
        str(k): int(v)
        for k, v in audit.get("root_p90_owner_counts", {}).items()
    }
    root_category_counts = {
        str(k): int(v)
        for k, v in audit.get("root_p90_category_counts", {}).items()
    }
    residual_items = tuple(str(item) for item in audit.get(
        "root_p90_current_residual_items", ()))
    blockers = tuple(str(item) for item in decision.get("promotion_blockers", ()))
    owner_counts_2500: dict[str, int] = {}
    for row in ci_worlds:
        if int(row.get("cells", 0)) == 2500:
            owner_counts_2500 = {
                str(k): int(v)
                for k, v in row.get("owner_counts", {}).items()
            }
            break
    packets = _packet_blueprints(
        owner_counts_900=root_owner_counts,
        owner_counts_2500=owner_counts_2500,
        category_counts=root_category_counts,
        residual_items=residual_items,
        blockers=blockers,
    )
    packet_ids = {packet["packet_id"] for packet in packets}
    blocker_assignments = {
        blocker: tuple(packet for packet in BLOCKER_TO_PACKET.get(blocker, ())
                       if packet in packet_ids)
        for blocker in blockers
    }
    residual_item_assignments = {
        item: OWNER_TO_PACKET.get(_owner_for_residual_item(item), "")
        for item in residual_items
    }
    unassigned_blockers = tuple(
        blocker for blocker, assigned in blocker_assignments.items()
        if not assigned
    )
    unassigned_residual_items = _unassigned_residual_items(residual_items)
    owner_layers = tuple(sorted(root_owner_counts))
    unassigned_owner_layers = tuple(
        owner for owner in owner_layers
        if owner not in OWNER_TO_PACKET
    )
    dependency_order_valid = _dependency_order_valid(packets)
    climate_targets = tuple(
        target
        for packet in packets
        for target in packet["implementation_targets"]
        if "climate" in str(target).lower()
        or "monsoon" in str(target).lower()
        or "ocean-current" in str(target).lower()
    )
    metrics = {
        "p91_summary_available": bool(p91_summary),
        "p91_status": str(p91_summary.get("status", "")),
        "p91_audit_status": str(audit.get("status", "")),
        "p91_promotion_ready": bool(decision.get("promotion_ready", False)),
        "p91_audit_completed": bool(decision.get("audit_completed", False)),
        "promotion_blocker_count": int(len(blockers)),
        "assigned_blocker_count": int(
            len(blockers) - len(unassigned_blockers)),
        "unassigned_blocker_count": int(len(unassigned_blockers)),
        "owner_layer_count": int(len(owner_layers)),
        "assigned_owner_layer_count": int(
            len(owner_layers) - len(unassigned_owner_layers)),
        "unassigned_owner_layer_count": int(len(unassigned_owner_layers)),
        "residual_item_count": int(len(residual_items)),
        "assigned_residual_item_count": int(
            len(residual_items) - len(unassigned_residual_items)),
        "unassigned_residual_item_count": int(len(unassigned_residual_items)),
        "repair_packet_count": int(len(packets)),
        "dependency_order_valid": bool(dependency_order_valid),
        "climate_target_count": int(len(climate_targets)),
        "packets_with_microbenchmarks": int(sum(
            bool(packet.get("microbenchmarks_to_add")) for packet in packets)),
        "packets_with_acceptance_targets": int(sum(
            bool(packet.get("acceptance_targets")) for packet in packets)),
        "packets_with_validation_suites": int(sum(
            bool(packet.get("validation_suites")) for packet in packets)),
        "packets_with_implementation_targets": int(sum(
            bool(packet.get("implementation_targets")) for packet in packets)),
        "final_reaudit_packet_defined": bool(
            packets and packets[-1]["packet_id"] == (
                "P92.8_integrated_reaudit_and_promotion_gate")),
        "next_implementation_packet": packets[0]["packet_id"] if packets else "",
        "final_validation_suite": "P91",
    }
    acceptance = {
        "p91_audit_available": bool(
            metrics["p91_summary_available"]
            and metrics["p91_status"] == "pass"
            and metrics["p91_audit_status"]
            == "integrated_real_earth_morphology_promotion_audit_ready"
        ),
        "p91_promotion_blocked_not_ignored": bool(
            not metrics["p91_promotion_ready"]
            and metrics["promotion_blocker_count"] > 0
        ),
        "all_p91_blockers_assigned": metrics["unassigned_blocker_count"] == 0,
        "all_p90_owner_layers_assigned": metrics["unassigned_owner_layer_count"] == 0,
        "all_residual_items_assigned": metrics["unassigned_residual_item_count"] == 0,
        "repair_packets_have_microbenchmarks": (
            metrics["packets_with_microbenchmarks"] == metrics["repair_packet_count"]),
        "repair_packets_have_acceptance_targets": (
            metrics["packets_with_acceptance_targets"] == metrics["repair_packet_count"]),
        "repair_packets_have_validation_suites": (
            metrics["packets_with_validation_suites"] == metrics["repair_packet_count"]),
        "repair_packets_have_implementation_targets": (
            metrics["packets_with_implementation_targets"] == metrics["repair_packet_count"]),
        "dependencies_are_ordered": metrics["dependency_order_valid"],
        "climate_ocean_monsoon_work_excluded": metrics["climate_target_count"] == 0,
        "final_p91_reaudit_defined": metrics["final_reaudit_packet_defined"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "production_residual_owner_repair_plan_ready"
            if all(acceptance.values())
            else "production_residual_owner_repair_plan_incomplete"
        ),
        "selected_p91_summary_path": selected_path,
        "promotion_blockers": blockers,
        "owner_layers": owner_layers,
        "root_p90_owner_counts": root_owner_counts,
        "root_p90_category_counts": root_category_counts,
        "root_p90_current_residual_items": residual_items,
        "repair_packets": tuple(packets),
        "blocker_assignments": blocker_assignments,
        "residual_item_assignments": residual_item_assignments,
        "unassigned_blockers": unassigned_blockers,
        "unassigned_owner_layers": unassigned_owner_layers,
        "unassigned_residual_items": unassigned_residual_items,
        "climate_targets": climate_targets,
        "metrics": metrics,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "packets": packets,
            "blockers": blockers,
            "items": residual_items,
            "metrics": metrics,
        }),
    }
