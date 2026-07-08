"""P93 planform reference calibration diagnostics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from aevum.diagnostics.earth_reference import REFERENCE_ENVELOPES


SCHEMA = "aevum.planform_reference_calibration.v1"

PLANFORM_CALIBRATION_METRICS = (
    "land_fraction",
    "largest_land_component_fraction",
    "land_component_count",
    "land_ribbon_fraction_gt_0_5",
    "land_coastline_complexity_largest",
)

CROSS_OWNER_METRICS = {
    "trench_fraction_of_ocean": "P92.7_bathymetry_margin_sequence",
}

P93_MICROBENCHMARKS = (
    "P93.planform_reference_calibration",
    "P93.generated_component_ribbon_envelope",
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _latest_summary(root: Path, pattern: str, schema: str) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob(pattern)):
        try:
            loaded = json.loads(path.read_text())
        except Exception:
            continue
        if loaded.get("schema") == schema:
            rows.append({"path": str(path), "summary": loaded})
    passing = [row for row in rows if row["summary"].get("status") == "pass"]
    return (passing or rows)[-1] if rows else {"path": "", "summary": {}}


def _metric_value_from_p90_gap(gaps: list[dict[str, Any]], metric: str) -> float | None:
    prefix = f"planform.{metric}_out_of_envelope"
    for gap in gaps:
        if gap.get("gap_id") == prefix and isinstance(gap.get("metric_value"), (int, float)):
            return float(gap["metric_value"])
    return None


def _severity(value: float | None, lo: float, hi: float) -> float:
    if value is None:
        return 0.0
    if value < lo:
        return float(lo / max(value, 1.0e-12))
    if value > hi:
        return float(value / max(hi, 1.0e-12))
    return 1.0


def _direction(value: float | None, lo: float, hi: float) -> str:
    if value is None:
        return "unknown"
    if value < lo:
        return "increase"
    if value > hi:
        return "decrease"
    return "hold"


def _p69_member_planform_rows(p69_summary: dict[str, Any]) -> tuple[dict[str, Any], ...]:
    metrics = p69_summary.get("benchmarks", [{}])[0].get("metrics", {})
    rows: list[dict[str, Any]] = []
    for row in metrics.get("per_member_metrics", ()):
        rows.append({
            "member_name": str(row.get("member_name", "")),
            "land_fraction": float(row.get("land_fraction", 0.0)),
            "largest_land_component_fraction": float(
                row.get("largest_land_component_fraction", 0.0)),
            "land_component_count": float(row.get("land_component_count", 0.0)),
            "land_ribbon_fraction_gt_0_5": float(row.get("land_ribbon_fraction", 0.0)),
            "continental_ribbon_fraction_gt_0_5": float(
                row.get("continental_ribbon_fraction", 0.0)),
            "land_coastline_complexity_largest": float(
                row.get("land_coastline_complexity_largest", 0.0)),
            "earth_reference_out_of_envelope": tuple(
                row.get("earth_reference_out_of_envelope", ())),
            "earth_reference_out_of_envelope_count": int(
                row.get("earth_reference_out_of_envelope_count", 0)),
        })
    return tuple(rows)


def _member_values(rows: tuple[dict[str, Any], ...], metric: str) -> tuple[float, ...]:
    return tuple(float(row.get(metric, 0.0)) for row in rows)


def _target_status(value: float | None, lo: float, hi: float) -> str:
    if value is None:
        return "missing_value"
    return "in_envelope" if lo <= value <= hi else "out_of_envelope"


def planform_reference_calibration_summary(root: Path) -> dict[str, Any]:
    p92 = _latest_summary(
        root,
        "out_bench_p92_production_residual_owner_repair_plan_*/tectonics_bench_summary.json",
        "aevum.tectonics_bench.p92.v1",
    )
    p91 = _latest_summary(
        root,
        "out_bench_p91_integrated_real_earth_morphology_promotion_audit_*/tectonics_bench_summary.json",
        "aevum.tectonics_bench.p91.v1",
    )
    p90 = _latest_summary(
        root,
        "out_bench_p90_current_world_morphology_gap_inventory_*/tectonics_bench_summary.json",
        "aevum.tectonics_bench.p90.v1",
    )
    p78 = _latest_summary(
        root,
        "out_bench_p78_generated_hypsometry_envelope_*/tectonics_bench_summary.json",
        "aevum.tectonics_bench.p78.v1",
    )
    p69 = _latest_summary(
        root,
        "out_bench_p69_highres_physical_ensemble_visual_audit_*/tectonics_bench_summary.json",
        "aevum.tectonics_bench.p69.v1",
    )

    p92_summary = p92["summary"]
    p91_summary = p91["summary"]
    p90_summary = p90["summary"]
    p78_summary = p78["summary"]
    p69_summary = p69["summary"]

    p92_plan = p92_summary.get("production_residual_owner_repair_plan", {})
    repair_packets = list(p92_plan.get("repair_packets", ()))
    first_packet = repair_packets[0] if repair_packets else {}
    p90_inventory = p90_summary.get("current_world_morphology_gap_inventory", {})
    p90_gaps = list(p90_inventory.get("gaps", ()))
    p90_planform_gaps = tuple(
        gap for gap in p90_gaps
        if gap.get("owner_layer") == "planform"
    )
    p78_metrics = p78_summary.get("benchmarks", [{}])[0].get("metrics", {})
    p69_rows = _p69_member_planform_rows(p69_summary)
    p69_earthlike = next(
        (row for row in p69_rows if row["member_name"] == "earthlike_reference_physical_member"),
        {},
    )
    p91_audit = p91_summary.get("integrated_real_earth_morphology_promotion_audit", {})
    p91_decision = p91_audit.get("promotion_decision", {})

    calibration_targets: list[dict[str, Any]] = []
    generated_envelope: dict[str, dict[str, Any]] = {}
    for metric in PLANFORM_CALIBRATION_METRICS:
        spec = REFERENCE_ENVELOPES[metric]
        lo, hi = (float(spec["range"][0]), float(spec["range"][1]))
        current_900 = _metric_value_from_p90_gap(p90_planform_gaps, metric)
        p69_values = _member_values(p69_rows, metric)
        earthlike_value = (
            float(p69_earthlike.get(metric, 0.0))
            if p69_earthlike else None
        )
        generated_envelope[metric] = {
            "reference_min": lo,
            "reference_max": hi,
            "current_900_value": current_900,
            "current_900_status": _target_status(current_900, lo, hi),
            "p69_member_min": min(p69_values) if p69_values else 0.0,
            "p69_member_max": max(p69_values) if p69_values else 0.0,
            "p69_member_in_envelope_count": int(sum(lo <= value <= hi for value in p69_values)),
            "p69_member_count": int(len(p69_values)),
            "p69_earthlike_value": earthlike_value,
            "p69_earthlike_status": _target_status(earthlike_value, lo, hi),
            "basis": spec["basis"],
        }
        calibration_targets.append({
            "metric": metric,
            "owner_layer": "planform",
            "reference_min": lo,
            "reference_max": hi,
            "current_900_value": current_900,
            "current_900_direction": _direction(current_900, lo, hi),
            "current_900_severity": _severity(current_900, lo, hi),
            "p69_earthlike_value": earthlike_value,
            "p69_earthlike_direction": _direction(earthlike_value, lo, hi),
            "p69_earthlike_severity": _severity(earthlike_value, lo, hi),
            "repair_packet": "P92.1_planform_and_reference_calibration",
            "production_knobs": _production_knobs_for_metric(metric),
            "target_microbenchmarks": P93_MICROBENCHMARKS,
        })

    cross_owner_targets: list[dict[str, Any]] = []
    for metric, packet in CROSS_OWNER_METRICS.items():
        if metric in REFERENCE_ENVELOPES:
            spec = REFERENCE_ENVELOPES[metric]
            lo, hi = (float(spec["range"][0]), float(spec["range"][1]))
        else:
            lo, hi = 0.0, 0.0
        current_900 = _metric_value_from_p90_gap(p90_planform_gaps, metric)
        cross_owner_targets.append({
            "metric": metric,
            "owner_layer": "bathymetry_margin",
            "reference_min": lo,
            "reference_max": hi,
            "current_900_value": current_900,
            "current_900_direction": _direction(current_900, lo, hi),
            "current_900_severity": _severity(current_900, lo, hi),
            "repair_packet": packet,
            "reason": "P90 records this as a planform envelope gap, but the production repair belongs to the margin/bathymetry sequence.",
        })

    p90_planform_gap_ids = tuple(str(gap.get("gap_id", "")) for gap in p90_planform_gaps)
    planform_metrics_covered = tuple(
        target["metric"] for target in calibration_targets
        if any(str(target["metric"]) in gap_id for gap_id in p90_planform_gap_ids)
        or target["metric"] == "largest_land_component_fraction"
    )
    unresolved_primary = tuple(
        target["metric"] for target in calibration_targets
        if target["current_900_direction"] != "hold"
        or target["p69_earthlike_direction"] != "hold"
    )
    microbenchmarks_declared = tuple(first_packet.get("microbenchmarks_to_add", ()))

    metrics = {
        "p92_packet_available": bool(first_packet),
        "p92_packet_id": str(first_packet.get("packet_id", "")),
        "p91_promotion_ready": bool(p91_decision.get("promotion_ready", False)),
        "p91_blocker_count": int(len(p91_decision.get("promotion_blockers", ()))),
        "p90_planform_gap_count": int(len(p90_planform_gaps)),
        "p90_planform_wrong_area_scale_count": int(sum(
            gap.get("category") == "wrong_area_scale" for gap in p90_planform_gaps)),
        "calibration_target_count": int(len(calibration_targets)),
        "cross_owner_target_count": int(len(cross_owner_targets)),
        "covered_planform_metric_count": int(len(planform_metrics_covered)),
        "p69_member_count": int(len(p69_rows)),
        "p69_earthlike_out_of_envelope_count": int(
            p69_earthlike.get("earth_reference_out_of_envelope_count", 0)
            if p69_earthlike else 0),
        "p78_current_out_of_envelope_count": int(
            p78_metrics.get("current_out_of_envelope_count", 0)),
        "p78_archived_highres_out_of_envelope_max": int(
            p78_metrics.get("archived_highres_out_of_envelope_max", 0)),
        "unresolved_primary_planform_metric_count": int(len(unresolved_primary)),
        "microbenchmark_count": int(len(microbenchmarks_declared)),
        "generated_envelope_metric_count": int(len(generated_envelope)),
        "planform_blocker_preserved": bool(
            "planform_residuals_unresolved" in p91_decision.get("promotion_blockers", ())),
        "p69_calibration_blocker_preserved": bool(
            "p69_earthlike_reference_needs_calibration"
            in p91_decision.get("promotion_blockers", ())),
        "does_not_promote_default": bool(not p91_decision.get("promotion_ready", True)),
        "next_packet_after_p93": "P92.2_production_province_graph_fields",
    }
    acceptance = {
        "p92_planform_packet_available": (
            metrics["p92_packet_id"] == "P92.1_planform_and_reference_calibration"),
        "p69_p78_p90_evidence_available": bool(
            p69_summary and p78_summary and p90_summary),
        "p91_blockers_preserved": bool(
            metrics["planform_blocker_preserved"]
            and metrics["p69_calibration_blocker_preserved"]
            and metrics["does_not_promote_default"]),
        "all_planform_reference_ranges_available": all(
            metric in REFERENCE_ENVELOPES for metric in PLANFORM_CALIBRATION_METRICS),
        "p90_planform_gaps_covered": (
            metrics["p90_planform_gap_count"] >= 5
            and metrics["covered_planform_metric_count"] >= 4
            and metrics["cross_owner_target_count"] >= 1),
        "calibration_directions_defined": all(
            target["current_900_direction"] != "unknown"
            or target["p69_earthlike_direction"] != "unknown"
            for target in calibration_targets),
        "generated_900_and_8000_envelope_defined": (
            metrics["p69_member_count"] >= 3
            and metrics["generated_envelope_metric_count"] == len(PLANFORM_CALIBRATION_METRICS)),
        "candidate_microbenchmarks_declared": (
            set(microbenchmarks_declared) == set(P93_MICROBENCHMARKS)),
        "trench_gap_deferred_to_bathymetry_packet": any(
            target["metric"] == "trench_fraction_of_ocean"
            and target["repair_packet"] == "P92.7_bathymetry_margin_sequence"
            for target in cross_owner_targets),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "planform_reference_calibration_ready"
            if all(acceptance.values())
            else "planform_reference_calibration_incomplete"
        ),
        "source_summaries": {
            "P69": p69["path"],
            "P78": p78["path"],
            "P90": p90["path"],
            "P91": p91["path"],
            "P92": p92["path"],
        },
        "p92_planform_packet": first_packet,
        "p90_planform_gap_ids": p90_planform_gap_ids,
        "p78_current_out_of_envelope": tuple(
            p78_metrics.get("current_out_of_envelope", ())),
        "p69_planform_members": p69_rows,
        "generated_envelope": generated_envelope,
        "calibration_targets": tuple(calibration_targets),
        "cross_owner_targets": tuple(cross_owner_targets),
        "unresolved_primary_planform_metrics": unresolved_primary,
        "microbenchmarks_declared": microbenchmarks_declared,
        "metrics": metrics,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "targets": calibration_targets,
            "cross": cross_owner_targets,
            "metrics": metrics,
        }),
    }


def _production_knobs_for_metric(metric: str) -> tuple[str, ...]:
    knobs = {
        "land_fraction": (
            "target_land_fraction calibration",
            "late-stage exposed land conservation",
            "coastal progradation / drowning balance",
        ),
        "largest_land_component_fraction": (
            "continent split/merge lifecycle",
            "supercontinent allowance by preset",
        ),
        "land_component_count": (
            "major component preservation",
            "sliver filtering without component collapse",
        ),
        "land_ribbon_fraction_gt_0_5": (
            "continental ribbon suppression",
            "broad interior support",
            "parented sliver policy",
        ),
        "land_coastline_complexity_largest": (
            "coastline smoothing guard",
            "valid bay/island-arc preservation",
        ),
    }
    return knobs.get(metric, ())
