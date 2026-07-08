"""P110B seed-sweep diagnostics for Earth-like terminal planforms.

P110B is a distribution-level review over already-generated P107 audits.  It
does not run the world generator; it reads P107 terminal metrics and asks
whether a seed sweep clusters near a modern Earth-like multipolar planform or
near the upper supercontinent-like warning band.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCHEMA = "aevum.p110b_seed_sweep_summary.v1"


@dataclass(frozen=True)
class P110BSeedSweepThresholds:
    min_sample_size: int = 5
    max_soft_warning_rate: float = 0.25
    max_median_largest_land_share: float = 0.58
    max_p90_largest_land_share: float = 0.62
    soft_largest_land_share: float = 0.60
    preferred_largest_land_share: float = 0.56
    earth_reference_largest_land_share: float = 0.5404
    min_major_land_components: int = 4
    max_major_land_components: int = 7
    min_second_land_share: float = 0.10
    min_third_land_share: float = 0.07
    min_major_ocean_basins: int = 2
    max_largest_ocean_basin_share: float = 0.65
    max_closed_ocean_ring_score: float = 0.03
    max_historical_supercontinent_frame_fraction: float = 0.45
    max_historical_supercontinent_duration_myr: float = 1200.0
    p111_preferred_max_largest_land_share: float = 0.45
    p111_soft_max_largest_land_share: float = 0.55
    p111_min_major_land_components: int = 3
    p111_max_major_land_components: int = 6
    p111_min_third_land_share: float = 0.10
    p111_min_fourth_land_share: float = 0.03
    p111_min_secondary_pair_share: float = 0.16
    p111_max_two_body_top2_share: float = 0.88


def summarize_p110b_seed_sweep(
    inputs: Iterable[str | Path],
    *,
    thresholds: P110BSeedSweepThresholds | None = None,
) -> dict[str, Any]:
    """Read P107 audit JSON artefacts and return a P110B distribution report."""
    thresholds = thresholds or P110BSeedSweepThresholds()
    input_paths = tuple(Path(path) for path in inputs)
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    for source in _discover_metric_sources(input_paths):
        try:
            rows.extend(_load_source_rows(source, thresholds))
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            warnings.append(f"{source}: {exc}")
    rows.sort(key=lambda row: (
        int(row.get("cells") or 0),
        int(row.get("n_plates") or 0),
        -1 if row.get("seed") is None else int(row["seed"]),
        str(row.get("source", "")),
    ))
    aggregate = _aggregate_rows(rows, thresholds)
    return {
        "schema": SCHEMA,
        "thresholds": _threshold_dict(thresholds),
        "input_paths": [str(path) for path in input_paths],
        "source_count": int(len(_discover_metric_sources(input_paths))),
        "run_count": int(len(rows)),
        "runs": rows,
        "aggregate": aggregate,
        "warnings": warnings,
        "acceptance": _acceptance(aggregate, thresholds),
    }


def write_p110b_seed_sweep_summary(
    summary: dict[str, Any],
    out: str | Path,
) -> Path:
    """Write a P110B summary to either a JSON path or output directory."""
    out_path = Path(out)
    if out_path.suffix.lower() == ".json":
        target = out_path
        target.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_path.mkdir(parents=True, exist_ok=True)
        target = out_path / "p110b_seed_sweep_summary.json"
    target.write_text(json.dumps(summary, indent=2, default=_json_default) + "\n")
    return target


def _discover_metric_sources(paths: Iterable[Path]) -> list[Path]:
    sources: list[Path] = []
    seen: set[Path] = set()
    for raw in paths:
        path = raw.expanduser()
        if path.is_file():
            candidates = [path]
        elif path.is_dir():
            summary = path / "p107_audit_summary.json"
            if summary.exists():
                candidates = [summary]
            else:
                candidates = sorted(path.glob("**/p107_terminal_metrics.json"))
        else:
            candidates = []
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved not in seen:
                sources.append(candidate)
                seen.add(resolved)
    return sources


def _load_source_rows(
    source: Path,
    thresholds: P110BSeedSweepThresholds,
) -> list[dict[str, Any]]:
    data = json.loads(source.read_text())
    if _looks_like_p107_ladder(data):
        return _load_ladder_rows(source, data, thresholds)
    if _looks_like_terminal_metrics(data):
        return [_normalize_run_row(
            metrics=data,
            metadata={
                "label": source.parent.name,
                "preset": None,
                "cells": data.get("cells"),
                "n_plates": None,
                "seed": None,
                "outdir": str(source.parent),
                "source": str(source),
            },
            thresholds=thresholds,
        )]
    raise ValueError("not a P107 audit summary or terminal metrics file")


def _load_ladder_rows(
    source: Path,
    data: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_dir = source.parent
    for entry in data.get("entries", []):
        raw_outdir = str(entry.get("outdir") or "")
        outdir = Path(raw_outdir) if raw_outdir else Path()
        if raw_outdir and not outdir.is_absolute() and not outdir.exists():
            outdir = base_dir / outdir
        metrics = entry.get("metrics", {})
        terminal: Path | None = (
            outdir / "p107_terminal_metrics.json" if raw_outdir else None)
        if terminal is not None and terminal.exists():
            metrics = json.loads(terminal.read_text())
        metadata = {
            "label": entry.get("label"),
            "preset": entry.get("preset"),
            "cells": entry.get("cells"),
            "n_plates": entry.get("n_plates"),
            "seed": entry.get("seed"),
            "outdir": str(outdir) if str(outdir) else entry.get("outdir", ""),
            "source": str(
                terminal
                if terminal is not None and terminal.exists()
                else source
            ),
            "run_seconds": entry.get("run_seconds"),
        }
        rows.append(_normalize_run_row(
            metrics=metrics,
            metadata=metadata,
            thresholds=thresholds,
        ))
    return rows


def _normalize_run_row(
    *,
    metrics: dict[str, Any],
    metadata: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> dict[str, Any]:
    planform = metrics.get("p110a_modern_planform", {})
    plan_summary = planform.get("summary", {})
    archetype = planform.get("p110b_final_state_archetype", {})
    lineage = planform.get("p110b_lineage_survival", {})
    gateway = planform.get("ocean_gateway_topology", {})
    supercontinent = planform.get("p110b_terminal_supercontinent_diagnostics", {})
    seaway = planform.get("seaway_cut_effectiveness", {})
    history = metrics.get("p110b_historical_supercontinent_trajectory", {})
    land_topology = planform.get("land", {})
    top_land = land_topology.get("top_component_shares_of_mask", []) or []
    fourth_share = float(top_land[3]) if len(top_land) > 3 else 0.0
    p109 = metrics.get("p109_hypsometry_comparison", {})
    ocean_basins = planform.get("ocean_basins", {})
    warnings = tuple(str(x) for x in planform.get("warning_flags", []) or [])
    out_of_envelope = tuple(str(x) for x in planform.get("out_of_envelope", []) or [])
    p109_out = tuple(str(x) for x in p109.get("out_of_envelope", []) or [])
    lineage_major_count = _int(lineage.get(
        "major_component_count",
        plan_summary.get("major_land_component_count"),
    ))
    lineage_supported_count = _int(
        lineage.get("lineage_supported_major_component_count"))
    # Legacy v1 lineage metrics treated every major exposed land component as a
    # continental-lineage subject.  Keep old reports comparable while v2 reports
    # can separate oceanic exposed landforms from continental lineage subjects.
    continental_major_count = _int(lineage.get(
        "continental_major_component_count",
        lineage_major_count,
    ))
    continental_supported_count = _int(lineage.get(
        "continental_lineage_supported_major_component_count",
        lineage_supported_count,
    ))
    oceanic_landform_count = _int(
        lineage.get("oceanic_landform_major_component_count", 0))
    provenance_supported_count = _int(lineage.get(
        "provenance_supported_major_component_count",
        lineage_supported_count,
    ))

    row = {
        "label": metadata.get("label"),
        "preset": metadata.get("preset"),
        "cells": _maybe_int(metadata.get("cells", metrics.get("cells"))),
        "n_plates": _maybe_int(metadata.get("n_plates")),
        "seed": _maybe_int(metadata.get("seed")),
        "run_seconds": _maybe_float(metadata.get("run_seconds")),
        "outdir": str(metadata.get("outdir") or ""),
        "source": str(metadata.get("source") or ""),
        "land_fraction": _float(plan_summary.get(
            "land_fraction", metrics.get("land_fraction"))),
        "p109_pass": bool(p109.get("within_p109_envelope", False)),
        "p109_out_of_envelope": list(p109_out),
        "p110a_pass": bool(planform.get(
            "within_p110a_modern_planform_envelope", False)),
        "p110a_out_of_envelope": list(out_of_envelope),
        "p110a_warning_flags": list(warnings),
        "p110b_archetype_code": _int(archetype.get("code")),
        "p110b_archetype": archetype.get("name", "unknown"),
        "p110b_archetype_largest_soft_ceiling": _float(
            archetype.get("largest_share_soft_ceiling")),
        "p110b_lineage_supported_major_component_count": lineage_supported_count,
        "p110b_continental_major_component_count": continental_major_count,
        "p110b_continental_lineage_supported_major_component_count": (
            continental_supported_count
        ),
        "p110b_oceanic_landform_major_component_count": oceanic_landform_count,
        "p110b_provenance_supported_major_component_count": provenance_supported_count,
        "p110b_independent_primary_lineage_count": _int(
            lineage.get("independent_primary_lineage_count")),
        "p110b_lineage_area_fraction_of_land": _float(
            lineage.get("lineage_area_fraction_of_land")),
        "p110b_geologic_support_fraction_of_land": _float(
            lineage.get("geologic_support_fraction_of_land")),
        "p110b_min_major_dominant_continent_share": _float(
            lineage.get("min_major_dominant_continent_share")),
        "p110b_min_major_continent_lineage_area_share": _float(
            lineage.get("min_major_continent_lineage_area_share")),
        "p110b_min_major_geologic_support_fraction": _float(
            lineage.get("min_major_geologic_support_fraction")),
        "p110b_terminal_supercontinent_score": _float(
            supercontinent.get(
                "terminal_supercontinent_score",
                plan_summary.get("terminal_supercontinent_score"),
            )
        ),
        "p110b_terminal_supercontinent_like": bool(
            supercontinent.get(
                "terminal_supercontinent_like",
                plan_summary.get("terminal_supercontinent_like", False),
            )
        ),
        "p110b_modern_multipolar_overconnected": bool(
            supercontinent.get("modern_multipolar_overconnected", False)
        ),
        "p110b_largest_land_significant_continent_domain_count": _int(
            supercontinent.get(
                "largest_land_significant_continent_domain_count",
                plan_summary.get(
                    "largest_land_significant_continent_domain_count"),
            )
        ),
        "p110b_largest_land_effective_continent_domain_count": _float(
            supercontinent.get(
                "largest_land_effective_continent_domain_count",
                plan_summary.get(
                    "largest_land_effective_continent_domain_count"),
            )
        ),
        "p110b_largest_land_robust_piece_count_after_neck_removal": _int(
            supercontinent.get(
                "largest_land_robust_piece_count_after_neck_removal",
                plan_summary.get(
                    "largest_land_robust_piece_count_after_neck_removal"),
            )
        ),
        "p110b_largest_land_bridge_candidate_fraction": _float(
            supercontinent.get(
                "largest_land_bridge_candidate_fraction",
                plan_summary.get("largest_land_bridge_candidate_fraction"),
            )
        ),
        "p110b_internal_domain_seaway_opening_count": _int(
            supercontinent.get(
                "internal_domain_seaway_opening_count",
                plan_summary.get("internal_domain_seaway_opening_count"),
            )
        ),
        "p110b_internal_domain_seaway_area_fraction_world": _float(
            supercontinent.get(
                "internal_domain_seaway_area_fraction_world",
                plan_summary.get("internal_domain_seaway_area_fraction_world"),
            )
        ),
        "p110b_internal_domain_boundary_count": _int(
            supercontinent.get(
                "internal_domain_boundary_count",
                plan_summary.get("internal_domain_boundary_count"),
            )
        ),
        "p110b_internal_domain_boundary_area_fraction_world": _float(
            supercontinent.get(
                "internal_domain_boundary_area_fraction_world",
                plan_summary.get("internal_domain_boundary_area_fraction_world"),
            )
        ),
        "p1114_modern_endpoint_seaway_count": _int(
            seaway.get("p1114_modern_endpoint_seaway_count", 0)),
        "p1114_modern_endpoint_seaway_area_fraction_world": _float(
            seaway.get("p1114_modern_endpoint_seaway_area_fraction_world", 0.0)),
        "p1114_modern_endpoint_seaway_domain_backed_count": _int(
            seaway.get("p1114_modern_endpoint_seaway_domain_backed_count", 0)),
        "p110b_historical_supercontinent_frame_fraction": _float(
            history.get("supercontinent_frame_fraction")),
        "p110b_historical_supercontinent_max_duration_myr": _float(
            history.get("max_consecutive_supercontinent_duration_myr")),
        "p110b_historical_supercontinent_time_window_myr": _float(
            history.get("supercontinent_time_window_myr")),
        "p110b_historical_max_largest_land_component_share": _float(
            history.get("max_largest_land_component_share")),
        "p110b_historical_usable_frame_count": _int(
            history.get("usable_frame_count")),
        "p110b_historical_long_lived_supercontinent_like": bool(
            history.get("long_lived_supercontinent_like", False)),
        "major_land_component_count": _int(plan_summary.get(
            "major_land_component_count")),
        "largest_land_component_share": _float(plan_summary.get(
            "largest_land_component_share")),
        "second_land_component_share": _float(plan_summary.get(
            "second_land_component_share")),
        "third_land_component_share": _float(plan_summary.get(
            "third_land_component_share")),
        "fourth_land_component_share": _float(fourth_share),
        "major_ocean_basin_count": _int(plan_summary.get(
            "major_ocean_basin_count")),
        "p110b_terminal_ocean_gateway_count": _int(
            gateway.get(
                "terminal_gateway_count",
                plan_summary.get("terminal_ocean_gateway_count"),
            )
        ),
        "p110b_terminal_interbasin_ocean_gateway_count": _int(
            gateway.get(
                "terminal_interbasin_gateway_count",
                plan_summary.get("terminal_interbasin_ocean_gateway_count"),
            )
        ),
        "p110b_terminal_phase_backed_ocean_gateway_count": _int(
            gateway.get(
                "terminal_phase_backed_gateway_count",
                plan_summary.get("terminal_phase_backed_ocean_gateway_count"),
            )
        ),
        "p110b_terminal_ocean_gateway_system_count": _int(
            gateway.get(
                "terminal_gateway_system_count",
                plan_summary.get(
                    "terminal_ocean_gateway_system_count",
                    gateway.get(
                        "terminal_gateway_count",
                        plan_summary.get("terminal_ocean_gateway_count"),
                    ),
                ),
            )
        ),
        "p110b_terminal_interbasin_ocean_gateway_system_count": _int(
            gateway.get(
                "terminal_interbasin_gateway_system_count",
                plan_summary.get(
                    "terminal_interbasin_ocean_gateway_system_count",
                    gateway.get(
                        "terminal_interbasin_gateway_count",
                        plan_summary.get("terminal_interbasin_ocean_gateway_count"),
                    ),
                ),
            )
        ),
        "p110b_terminal_phase_backed_ocean_gateway_system_count": _int(
            gateway.get(
                "terminal_phase_backed_gateway_system_count",
                plan_summary.get(
                    "terminal_phase_backed_ocean_gateway_system_count",
                    gateway.get(
                        "terminal_phase_backed_gateway_count",
                        plan_summary.get("terminal_phase_backed_ocean_gateway_count"),
                    ),
                ),
            )
        ),
        "p110b_ocean_gateway_fragment_to_system_ratio": _float(
            gateway.get(
                "gateway_fragment_to_system_ratio",
                plan_summary.get("ocean_gateway_fragment_to_system_ratio", 1.0),
            )
        ),
        "p110b_tectonic_ocean_gateway_count": _int(
            gateway.get(
                "tectonic_gateway_count",
                plan_summary.get("tectonic_ocean_gateway_count"),
            )
        ),
        "p110b_restricted_ocean_fraction": _float(
            gateway.get(
                "restricted_ocean_fraction",
                plan_summary.get("restricted_ocean_fraction"),
            )
        ),
        "p110b_unbacked_major_disconnected_ocean_component_count": _int(
            gateway.get(
                "unbacked_major_disconnected_ocean_component_count",
                plan_summary.get(
                    "unbacked_major_disconnected_ocean_component_count"),
            )
        ),
        "largest_ocean_basin_share": _float(plan_summary.get(
            "largest_ocean_basin_share")),
        "closed_ocean_ring_score": _float(plan_summary.get(
            "closed_ocean_ring_score")),
        "ocean_basin_source": ocean_basins.get("source"),
        "object_backed_ocean_basins": ocean_basins.get("source") == "ocean_basin_objects",
        "island_arc_chain_count": _int(metrics.get("island_arc_chain_count")),
        "microcontinent_object_count": _int(metrics.get("microcontinent_object_count")),
        "parented_oceanic_island_chain_count": _int(
            metrics.get("parented_oceanic_island_chain_count")),
        "deep_trench_fraction_below_6000m": _float(
            metrics.get("deep_trench_fraction_below_6000m")),
    }
    if _modern_balanced_three_major(row, thresholds):
        p110a_out = set(row["p110a_out_of_envelope"])
        if not row["p110a_pass"] and p110a_out <= {"major_land_component_count_low"}:
            row["p110a_pass"] = True
            row["p110a_out_of_envelope"] = []
            row["p110a_reclassified_modern_balanced_three_major"] = True
            if (
                "major_land_component_count_low_p110b_modern_balanced_three_major"
                not in row["p110a_warning_flags"]
            ):
                row["p110a_warning_flags"].append(
                    "major_land_component_count_low_p110b_modern_balanced_three_major"
                )
        else:
            row["p110a_reclassified_modern_balanced_three_major"] = False
    else:
        row["p110a_reclassified_modern_balanced_three_major"] = False
    status, review_flags = _classify_row(row, thresholds)
    row["p110b_status"] = status
    row["p110b_review_flags"] = review_flags
    p111_flags = _p111_modern_planform_flags(row, thresholds)
    row["p111_modern_planform_status"] = (
        "p111_modern_planform_candidate"
        if status == "p110b_visual_candidate" and not p111_flags
        else "p111_modern_planform_review"
    )
    row["p111_modern_planform_review_flags"] = p111_flags
    return row


def _classify_row(
    row: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> tuple[str, list[str]]:
    flags: list[str] = []
    if not row["p109_pass"]:
        flags.append("p109_hypsometry_failure")
    if not row["p110a_pass"]:
        flags.append("p110a_planform_failure")
    post_breakup_three_major = (
        row["p110b_archetype"] == "post_supercontinent_breakup"
        and row["major_land_component_count"] >= 3
        and row["second_land_component_share"] >= 0.16
        and row["third_land_component_share"] >= 0.16
    )
    archipelago_two_major = (
        row["p110b_archetype"] == "archipelago_active_margin"
        and row["major_land_component_count"] >= 2
        and row["second_land_component_share"] >= 0.20
        and row["island_arc_chain_count"] >= 12
        and row["microcontinent_object_count"] >= 6
    )
    modern_balanced_three_major = _modern_balanced_three_major(row, thresholds)
    split_secondary_pair = (
        row["major_land_component_count"] >= 4
        and row["second_land_component_share"] >= 0.18
        and row["third_land_component_share"] >= 0.040
        and row.get("fourth_land_component_share", 0.0) >= 0.030
        and (
            row["third_land_component_share"]
            + row.get("fourth_land_component_share", 0.0)
        ) >= 0.080
    )
    if row["major_land_component_count"] < thresholds.min_major_land_components:
        if not (
            post_breakup_three_major
            or archipelago_two_major
            or modern_balanced_three_major
        ):
            flags.append("major_land_component_count_low")
    if row["major_land_component_count"] > thresholds.max_major_land_components:
        flags.append("major_land_component_count_high")
    if row["largest_land_component_share"] > thresholds.soft_largest_land_share:
        flags.append("largest_land_component_share_soft")
    elif row["largest_land_component_share"] > thresholds.preferred_largest_land_share:
        flags.append("largest_land_component_share_above_preferred")
    if row["second_land_component_share"] < thresholds.min_second_land_share:
        flags.append("second_land_component_share_low")
    if (
        row["third_land_component_share"] < thresholds.min_third_land_share
        and not archipelago_two_major
        and not split_secondary_pair
    ):
        flags.append("third_land_component_share_low")
    if row["major_ocean_basin_count"] < thresholds.min_major_ocean_basins:
        flags.append("major_ocean_basin_count_low")
    if row["largest_ocean_basin_share"] > thresholds.max_largest_ocean_basin_share:
        flags.append("largest_ocean_basin_share_high")
    if row["closed_ocean_ring_score"] > thresholds.max_closed_ocean_ring_score:
        flags.append("closed_ocean_ring_score_high")
    if row["p110b_terminal_supercontinent_like"]:
        flags.append("terminal_supercontinent_like")
    if row["p110b_modern_multipolar_overconnected"]:
        flags.append("modern_multipolar_overconnected")
    if row["p110b_historical_long_lived_supercontinent_like"]:
        flags.append("historical_long_lived_supercontinent_like")
    if row["ocean_basin_source"] is None:
        flags.append("ocean_basin_source_missing")
    elif not row["object_backed_ocean_basins"]:
        flags.append("ocean_basin_source_not_object_backed")

    if not row["p109_pass"] or not row["p110a_pass"]:
        return "hard_gate_failure", flags
    if flags:
        return "threshold_pass_needs_p110b_review", flags
    return "p110b_visual_candidate", flags


def _modern_balanced_three_major(
    row: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> bool:
    """Return true for a modern, Earth-like three-major-landmass terminal state."""
    return bool(
        row["p110b_archetype"] == "modern_multipolar"
        and row["major_land_component_count"] >= 3
        and row["largest_land_component_share"] <= thresholds.earth_reference_largest_land_share
        and row["second_land_component_share"] >= 0.20
        and row["third_land_component_share"] >= 0.16
        and not row["p110b_terminal_supercontinent_like"]
        and not row["p110b_modern_multipolar_overconnected"]
        and (
            row["p110b_largest_land_significant_continent_domain_count"] >= 2
            or row["p110b_independent_primary_lineage_count"] >= 2
        )
    )


def _p111_modern_planform_flags(
    row: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> list[str]:
    """Return stricter P111 modern-endpoint review flags.

    P110B answers whether a generated map avoids the coarse supercontinent
    failure band.  P111 is stricter: it asks whether the terminal frame has a
    modern Earth-like distribution of several meaningful landmasses and open
    ocean basins.  Keep these flags separate so historical P110B reports remain
    comparable.
    """
    flags: list[str] = []
    major_count = int(row["major_land_component_count"])
    largest = float(row["largest_land_component_share"])
    second = float(row["second_land_component_share"])
    third = float(row["third_land_component_share"])
    fourth = float(row.get("fourth_land_component_share", 0.0))
    top2 = largest + second

    if major_count < thresholds.p111_min_major_land_components:
        flags.append("p111_major_land_component_count_low")
    if (
        major_count <= 2
        or (top2 >= thresholds.p111_max_two_body_top2_share and third < 0.07)
    ):
        flags.append("p111_two_body_endpoint")
    if major_count > thresholds.p111_max_major_land_components:
        flags.append("p111_major_land_component_count_high")

    if largest > thresholds.p111_soft_max_largest_land_share:
        flags.append("p111_largest_land_component_overlarge")
    elif largest > thresholds.p111_preferred_max_largest_land_share:
        flags.append("p111_largest_land_component_above_modern_preferred")

    secondary_tail = third + fourth
    if (
        major_count >= thresholds.p111_min_major_land_components
        and third < thresholds.p111_min_third_land_share
        and (
            fourth < thresholds.p111_min_fourth_land_share
            or secondary_tail < thresholds.p111_min_secondary_pair_share
        )
    ):
        flags.append("p111_weak_third_fourth_land_component")

    supported_secondary_count = min(
        int(row["p110b_provenance_supported_major_component_count"]),
        int(row["p110b_continental_lineage_supported_major_component_count"]),
    )
    if major_count >= 3 and supported_secondary_count < 3:
        flags.append("p111_secondary_component_geologic_support_weak")

    if (
        largest > thresholds.p111_soft_max_largest_land_share
        and int(row["p110b_internal_domain_seaway_opening_count"]) <= 0
        and int(row.get("p1114_modern_endpoint_seaway_count", 0)) <= 0
        and int(row["p110b_largest_land_robust_piece_count_after_neck_removal"]) <= 1
    ):
        flags.append("p111_overlarge_landmass_lacks_internal_seaway_provenance")

    return flags


def _aggregate_rows(
    rows: list[dict[str, Any]],
    thresholds: P110BSeedSweepThresholds,
) -> dict[str, Any]:
    run_count = len(rows)
    p109_pass_count = sum(1 for row in rows if row["p109_pass"])
    p110a_pass_count = sum(1 for row in rows if row["p110a_pass"])
    threshold_pass_count = sum(
        1 for row in rows if row["p109_pass"] and row["p110a_pass"])
    visual_candidate_count = sum(
        1 for row in rows if row["p110b_status"] == "p110b_visual_candidate")
    soft_warning_count = sum(
        1 for row in rows
        if "largest_land_component_share_soft" in row["p110b_review_flags"]
    )
    above_preferred_count = sum(
        1 for row in rows
        if "largest_land_component_share_above_preferred" in row["p110b_review_flags"]
    )
    object_backed_count = sum(1 for row in rows if row["object_backed_ocean_basins"])
    p111_candidate_count = sum(
        1 for row in rows
        if row["p111_modern_planform_status"] == "p111_modern_planform_candidate"
    )
    p111_review_count = run_count - p111_candidate_count
    p111_two_body_count = sum(
        1 for row in rows
        if "p111_two_body_endpoint" in row["p111_modern_planform_review_flags"]
    )
    p111_overlarge_count = sum(
        1 for row in rows
        if (
            "p111_largest_land_component_overlarge"
            in row["p111_modern_planform_review_flags"]
        )
    )
    p111_above_preferred_count = sum(
        1 for row in rows
        if (
            "p111_largest_land_component_above_modern_preferred"
            in row["p111_modern_planform_review_flags"]
        )
    )
    p111_weak_third_fourth_count = sum(
        1 for row in rows
        if (
            "p111_weak_third_fourth_land_component"
            in row["p111_modern_planform_review_flags"]
        )
    )
    p111_seaway_gap_count = sum(
        1 for row in rows
        if (
            "p111_overlarge_landmass_lacks_internal_seaway_provenance"
            in row["p111_modern_planform_review_flags"]
        )
    )
    aggregate = {
        "p109_pass_count": int(p109_pass_count),
        "p110a_pass_count": int(p110a_pass_count),
        "threshold_pass_count": int(threshold_pass_count),
        "p110b_visual_candidate_count": int(visual_candidate_count),
        "largest_land_component_soft_warning_count": int(soft_warning_count),
        "largest_land_component_above_preferred_count": int(above_preferred_count),
        "object_backed_ocean_basin_count": int(object_backed_count),
        "rates": {
            "p109_pass_rate": _rate(p109_pass_count, run_count),
            "p110a_pass_rate": _rate(p110a_pass_count, run_count),
            "threshold_pass_rate": _rate(threshold_pass_count, run_count),
            "p110b_visual_candidate_rate": _rate(visual_candidate_count, run_count),
            "largest_land_component_soft_warning_rate": _rate(
                soft_warning_count, run_count),
            "object_backed_ocean_basin_rate": _rate(object_backed_count, run_count),
            "p111_modern_planform_candidate_rate": _rate(
                p111_candidate_count, run_count),
            "p111_modern_planform_review_rate": _rate(
                p111_review_count, run_count),
            "p111_two_body_endpoint_rate": _rate(p111_two_body_count, run_count),
            "p111_largest_land_component_overlarge_rate": _rate(
                p111_overlarge_count, run_count),
            "p111_largest_land_component_above_modern_preferred_rate": _rate(
                p111_above_preferred_count, run_count),
            "p111_weak_third_fourth_land_component_rate": _rate(
                p111_weak_third_fourth_count, run_count),
            "p111_internal_seaway_provenance_gap_rate": _rate(
                p111_seaway_gap_count, run_count),
        },
        "largest_land_component_share": _series_stats(
            row["largest_land_component_share"] for row in rows),
        "second_land_component_share": _series_stats(
            row["second_land_component_share"] for row in rows),
        "third_land_component_share": _series_stats(
            row["third_land_component_share"] for row in rows),
        "land_fraction": _series_stats(row["land_fraction"] for row in rows),
        "major_land_component_count": _series_stats(
            row["major_land_component_count"] for row in rows),
        "major_ocean_basin_count": _series_stats(
            row["major_ocean_basin_count"] for row in rows),
        "p110b_terminal_ocean_gateway_count": _series_stats(
            row["p110b_terminal_ocean_gateway_count"] for row in rows),
        "p110b_terminal_interbasin_ocean_gateway_count": _series_stats(
            row["p110b_terminal_interbasin_ocean_gateway_count"] for row in rows),
        "p110b_terminal_phase_backed_ocean_gateway_count": _series_stats(
            row["p110b_terminal_phase_backed_ocean_gateway_count"] for row in rows),
        "p110b_terminal_ocean_gateway_system_count": _series_stats(
            row["p110b_terminal_ocean_gateway_system_count"] for row in rows),
        "p110b_terminal_interbasin_ocean_gateway_system_count": _series_stats(
            row["p110b_terminal_interbasin_ocean_gateway_system_count"]
            for row in rows),
        "p110b_terminal_phase_backed_ocean_gateway_system_count": _series_stats(
            row["p110b_terminal_phase_backed_ocean_gateway_system_count"]
            for row in rows),
        "p110b_ocean_gateway_fragment_to_system_ratio": _series_stats(
            row["p110b_ocean_gateway_fragment_to_system_ratio"] for row in rows),
        "p110b_tectonic_ocean_gateway_count": _series_stats(
            row["p110b_tectonic_ocean_gateway_count"] for row in rows),
        "p110b_restricted_ocean_fraction": _series_stats(
            row["p110b_restricted_ocean_fraction"] for row in rows),
        "p110b_unbacked_major_disconnected_ocean_component_count": _series_stats(
            row["p110b_unbacked_major_disconnected_ocean_component_count"]
            for row in rows),
        "largest_ocean_basin_share": _series_stats(
            row["largest_ocean_basin_share"] for row in rows),
        "closed_ocean_ring_score": _series_stats(
            row["closed_ocean_ring_score"] for row in rows),
        "island_arc_chain_count": _series_stats(
            row["island_arc_chain_count"] for row in rows),
        "microcontinent_object_count": _series_stats(
            row["microcontinent_object_count"] for row in rows),
        "parented_oceanic_island_chain_count": _series_stats(
            row["parented_oceanic_island_chain_count"] for row in rows),
        "deep_trench_fraction_below_6000m": _series_stats(
            row["deep_trench_fraction_below_6000m"] for row in rows),
        "p110b_lineage_supported_major_component_count": _series_stats(
            row["p110b_lineage_supported_major_component_count"] for row in rows),
        "p110b_continental_major_component_count": _series_stats(
            row["p110b_continental_major_component_count"] for row in rows),
        "p110b_continental_lineage_supported_major_component_count": _series_stats(
            row["p110b_continental_lineage_supported_major_component_count"]
            for row in rows),
        "p110b_oceanic_landform_major_component_count": _series_stats(
            row["p110b_oceanic_landform_major_component_count"] for row in rows),
        "p110b_provenance_supported_major_component_count": _series_stats(
            row["p110b_provenance_supported_major_component_count"] for row in rows),
        "p110b_independent_primary_lineage_count": _series_stats(
            row["p110b_independent_primary_lineage_count"] for row in rows),
        "p110b_lineage_area_fraction_of_land": _series_stats(
            row["p110b_lineage_area_fraction_of_land"] for row in rows),
        "p110b_geologic_support_fraction_of_land": _series_stats(
            row["p110b_geologic_support_fraction_of_land"] for row in rows),
        "p110b_min_major_dominant_continent_share": _series_stats(
            row["p110b_min_major_dominant_continent_share"] for row in rows),
        "p110b_min_major_continent_lineage_area_share": _series_stats(
            row["p110b_min_major_continent_lineage_area_share"] for row in rows),
        "p110b_min_major_geologic_support_fraction": _series_stats(
            row["p110b_min_major_geologic_support_fraction"] for row in rows),
        "p110b_terminal_supercontinent_score": _series_stats(
            row["p110b_terminal_supercontinent_score"] for row in rows),
        "p110b_terminal_supercontinent_like_count": int(sum(
            1 for row in rows if row["p110b_terminal_supercontinent_like"])),
        "p110b_modern_multipolar_overconnected_count": int(sum(
            1 for row in rows if row["p110b_modern_multipolar_overconnected"])),
        "p110b_largest_land_significant_continent_domain_count": _series_stats(
            row["p110b_largest_land_significant_continent_domain_count"]
            for row in rows),
        "p110b_largest_land_effective_continent_domain_count": _series_stats(
            row["p110b_largest_land_effective_continent_domain_count"]
            for row in rows),
        "p110b_largest_land_robust_piece_count_after_neck_removal": _series_stats(
            row["p110b_largest_land_robust_piece_count_after_neck_removal"]
            for row in rows),
        "p110b_largest_land_bridge_candidate_fraction": _series_stats(
            row["p110b_largest_land_bridge_candidate_fraction"] for row in rows),
        "p110b_internal_domain_seaway_opening_count": _series_stats(
            row["p110b_internal_domain_seaway_opening_count"] for row in rows),
        "p110b_internal_domain_seaway_area_fraction_world": _series_stats(
            row["p110b_internal_domain_seaway_area_fraction_world"] for row in rows),
        "p110b_internal_domain_boundary_count": _series_stats(
            row["p110b_internal_domain_boundary_count"] for row in rows),
        "p110b_internal_domain_boundary_area_fraction_world": _series_stats(
            row["p110b_internal_domain_boundary_area_fraction_world"]
            for row in rows),
        "p1114_modern_endpoint_seaway_count": _series_stats(
            row["p1114_modern_endpoint_seaway_count"] for row in rows),
        "p1114_modern_endpoint_seaway_area_fraction_world": _series_stats(
            row["p1114_modern_endpoint_seaway_area_fraction_world"]
            for row in rows),
        "p1114_modern_endpoint_seaway_domain_backed_count": _series_stats(
            row["p1114_modern_endpoint_seaway_domain_backed_count"]
            for row in rows),
        "p110b_historical_supercontinent_frame_fraction": _series_stats(
            row["p110b_historical_supercontinent_frame_fraction"] for row in rows),
        "p110b_historical_supercontinent_max_duration_myr": _series_stats(
            row["p110b_historical_supercontinent_max_duration_myr"] for row in rows),
        "p110b_historical_supercontinent_time_window_myr": _series_stats(
            row["p110b_historical_supercontinent_time_window_myr"] for row in rows),
        "p110b_historical_max_largest_land_component_share": _series_stats(
            row["p110b_historical_max_largest_land_component_share"] for row in rows),
        "p110b_historical_usable_frame_count": _series_stats(
            row["p110b_historical_usable_frame_count"] for row in rows),
        "p110b_historical_long_lived_supercontinent_like_count": int(sum(
            1 for row in rows
            if row["p110b_historical_long_lived_supercontinent_like"])),
        "review_flag_counts": _flag_counts(rows),
        "p111_modern_planform_candidate_count": int(p111_candidate_count),
        "p111_modern_planform_review_count": int(p111_review_count),
        "p111_two_body_endpoint_count": int(p111_two_body_count),
        "p111_largest_land_component_overlarge_count": int(p111_overlarge_count),
        "p111_largest_land_component_above_modern_preferred_count": int(
            p111_above_preferred_count),
        "p111_weak_third_fourth_land_component_count": int(
            p111_weak_third_fourth_count),
        "p111_internal_seaway_provenance_gap_count": int(p111_seaway_gap_count),
        "p111_review_flag_counts": _p111_flag_counts(rows),
        "status_counts": _value_counts(row["p110b_status"] for row in rows),
        "p111_status_counts": _value_counts(
            row["p111_modern_planform_status"] for row in rows),
        "p110b_archetype_counts": _value_counts(
            row["p110b_archetype"] for row in rows),
        "major_land_component_count_histogram": _value_counts(
            str(row["major_land_component_count"]) for row in rows),
        "ocean_basin_source_counts": _value_counts(
            str(row["ocean_basin_source"]) for row in rows),
    }
    aggregate["distribution_flags"] = _distribution_flags(
        aggregate, run_count, thresholds)
    aggregate["p111_distribution_flags"] = _p111_distribution_flags(
        aggregate, run_count, thresholds)
    return aggregate


def _distribution_flags(
    aggregate: dict[str, Any],
    run_count: int,
    thresholds: P110BSeedSweepThresholds,
) -> list[str]:
    flags: list[str] = []
    rates = aggregate["rates"]
    largest = aggregate["largest_land_component_share"]
    if run_count < thresholds.min_sample_size:
        flags.append("sample_size_low")
    if aggregate["p109_pass_count"] != run_count:
        flags.append("p109_failures_present")
    if aggregate["p110a_pass_count"] != run_count:
        flags.append("p110a_failures_present")
    if rates["largest_land_component_soft_warning_rate"] > thresholds.max_soft_warning_rate:
        flags.append("largest_landmass_soft_warning_rate_high")
    if _stat_value(largest, "median") > thresholds.max_median_largest_land_share:
        flags.append("median_largest_landmass_too_high")
    if _stat_value(largest, "p90") > thresholds.max_p90_largest_land_share:
        flags.append("p90_largest_landmass_too_high")
    if aggregate["object_backed_ocean_basin_count"] != run_count:
        flags.append("object_backed_ocean_basin_source_incomplete")
    if aggregate.get("p110b_terminal_supercontinent_like_count", 0) > 0:
        flags.append("terminal_supercontinent_like_present")
    if aggregate.get("p110b_modern_multipolar_overconnected_count", 0) > 0:
        flags.append("modern_multipolar_overconnected_present")
    if aggregate.get("p110b_historical_long_lived_supercontinent_like_count", 0) > 0:
        flags.append("historical_long_lived_supercontinent_like_present")
    if aggregate["p110b_visual_candidate_count"] != run_count:
        flags.append("p110b_review_flags_present")
    return flags


def _p111_distribution_flags(
    aggregate: dict[str, Any],
    run_count: int,
    thresholds: P110BSeedSweepThresholds,
) -> list[str]:
    flags: list[str] = []
    largest = aggregate["largest_land_component_share"]
    major = aggregate["major_land_component_count"]
    rates = aggregate["rates"]
    if run_count < thresholds.min_sample_size:
        flags.append("sample_size_low")
    if aggregate["p111_modern_planform_candidate_count"] != run_count:
        flags.append("p111_modern_planform_review_flags_present")
    if aggregate["p111_two_body_endpoint_count"] > 0:
        flags.append("p111_two_body_endpoint_present")
    if aggregate["p111_largest_land_component_overlarge_count"] > 0:
        flags.append("p111_overlarge_largest_landmass_present")
    if aggregate["p111_weak_third_fourth_land_component_count"] > 0:
        flags.append("p111_weak_third_fourth_land_component_present")
    if aggregate["p111_internal_seaway_provenance_gap_count"] > 0:
        flags.append("p111_internal_seaway_provenance_gap_present")
    if _stat_value(largest, "median") > thresholds.p111_preferred_max_largest_land_share:
        flags.append("p111_median_largest_landmass_above_modern_preferred")
    if _stat_value(largest, "p90") > thresholds.p111_soft_max_largest_land_share:
        flags.append("p111_p90_largest_landmass_above_modern_soft")
    if _stat_value(major, "min") < thresholds.p111_min_major_land_components:
        flags.append("p111_major_land_component_count_low_present")
    if _stat_value(major, "max") > thresholds.p111_max_major_land_components:
        flags.append("p111_major_land_component_count_high_present")
    if rates.get("p111_largest_land_component_above_modern_preferred_rate", 0.0) > 0.5:
        flags.append("p111_largest_landmass_above_preferred_rate_high")
    return flags


def _acceptance(
    aggregate: dict[str, Any],
    thresholds: P110BSeedSweepThresholds,
) -> dict[str, Any]:
    run_count = int(sum(aggregate["status_counts"].values()))
    flags = set(aggregate["distribution_flags"])
    return {
        "p110b_seed_sweep_diagnostics_complete": True,
        "sample_size_sufficient": run_count >= thresholds.min_sample_size,
        "all_runs_pass_p109_and_p110a": (
            aggregate["p109_pass_count"] == run_count
            and aggregate["p110a_pass_count"] == run_count
        ),
        "object_basin_source_complete": (
            aggregate["object_backed_ocean_basin_count"] == run_count
        ),
        "soft_warning_rate_ok": (
            aggregate["rates"]["largest_land_component_soft_warning_rate"]
            <= thresholds.max_soft_warning_rate
        ),
        "median_largest_component_near_reference": (
            _stat_value(aggregate["largest_land_component_share"], "median")
            <= thresholds.max_median_largest_land_share
        ),
        "p90_largest_component_below_p110b_soft_ceiling": (
            _stat_value(aggregate["largest_land_component_share"], "p90")
            <= thresholds.max_p90_largest_land_share
        ),
        "historical_supercontinent_duration_ok": (
            "historical_long_lived_supercontinent_like_present" not in flags
        ),
        "ready_for_p110b_archetype_tuning": (
            run_count > 0 and "sample_size_low" not in flags
        ),
        "ready_for_p110b_promotion": len(flags) == 0,
        "p111_modern_planform_distribution_ready": (
            len(aggregate.get("p111_distribution_flags", [])) == 0
        ),
    }


def _looks_like_p107_ladder(data: dict[str, Any]) -> bool:
    return isinstance(data.get("entries"), list) and "p107_audit_ladder" in str(
        data.get("schema", ""))


def _looks_like_terminal_metrics(data: dict[str, Any]) -> bool:
    return isinstance(data.get("p110a_modern_planform"), dict)


def _series_stats(values: Iterable[Any]) -> dict[str, Any]:
    arr = np.asarray([float(x) for x in values if x is not None], dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "count": 0,
            "min": None,
            "median": None,
            "mean": None,
            "p90": None,
            "max": None,
        }
    return {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "p90": float(np.percentile(arr, 90.0)),
        "max": float(np.max(arr)),
    }


def _flag_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        for flag in row["p110b_review_flags"]:
            counts[flag] = counts.get(flag, 0) + 1
    return dict(sorted(counts.items()))


def _p111_flag_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        for flag in row["p111_modern_planform_review_flags"]:
            counts[flag] = counts.get(flag, 0) + 1
    return dict(sorted(counts.items()))


def _value_counts(values: Iterable[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _threshold_dict(thresholds: P110BSeedSweepThresholds) -> dict[str, Any]:
    return {
        key: getattr(thresholds, key)
        for key in thresholds.__dataclass_fields__
    }


def _rate(count: int, total: int) -> float:
    return float(count / total) if total else 0.0


def _stat_value(stats: dict[str, Any], key: str, default: float = 0.0) -> float:
    value = stats.get(key, default)
    if value is None:
        return float(default)
    return _float(value, default)


def _float(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    if not np.isfinite(out):
        return float(default)
    return out


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _maybe_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _maybe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)
