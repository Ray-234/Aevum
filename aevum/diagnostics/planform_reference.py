"""P30 planform diagnostics for Earth-like land and coastline geometry.

This module is a screening layer, not a direct ETOPO/GPlates fit.  It extracts
the current morphology metrics that explain the largest residual visual
problems: ribbon-like land and over-complex coastlines.
"""
from __future__ import annotations

from typing import Any

from aevum.diagnostics.earth_reference import (
    REFERENCE_ENVELOPES,
    REFERENCE_SOURCES,
)


SCHEMA = "aevum.planform_reference.v1"

PLANFORM_KEYS = (
    "land_fraction",
    "largest_land_component_fraction",
    "land_component_count",
    "land_ribbon_fraction_gt_0_5",
    "land_coastline_complexity_largest",
)

RESOLUTION_LADDER = [
    {
        "cells": "8000",
        "role": "routine Earth-like planform and topology review",
    },
    {
        "cells": "24000",
        "role": "medium audit for coastline, arc, shelf, strait, and delta geometry",
    },
    {
        "cells": "72000",
        "role": "occasional offline audit for isthmuses, compact mountains, narrow shelves, small deltas/fans, and short straits",
    },
]


def planform_reference_summary(
    world: Any | None = None,
    morphology_metrics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return Earth-reference screening metrics for continental planform.

    The summary intentionally uses broad envelopes from
    :mod:`aevum.diagnostics.earth_reference`.  It is meant to guide parameter
    tuning and regression gates before a direct raster/plate-model comparison is
    integrated.
    """
    if morphology_metrics is None:
        if world is None:
            raise ValueError("world or morphology_metrics must be supplied")
        from aevum.diagnostics.morphology import compute_world_morphology

        morphology_metrics = compute_world_morphology(world).metrics

    entry_land_fraction = None
    if "morphology" in morphology_metrics:
        entry_land_fraction = morphology_metrics.get("land_fraction")
        morphology_metrics = morphology_metrics.get("morphology", {})

    if "exposed_land" in morphology_metrics:
        exposed = morphology_metrics.get("exposed_land", {})
        continental = morphology_metrics.get("continental_crust", {})
        coupling = morphology_metrics.get("crust_land_coupling", {})
        values = {
            "land_fraction": float(
                coupling.get(
                    "land_area_fraction",
                    exposed.get("mask_area_fraction_of_total", 0.0),
                )
            ),
            "largest_land_component_fraction": float(
                exposed.get("largest_component_area_fraction_of_mask", 0.0)
            ),
            "land_component_count": float(exposed.get("component_count", 0.0)),
            "land_ribbon_fraction_gt_0_5": float(
                exposed.get("ribbon_area_fraction_gt_0_5", 0.0)
            ),
            "land_coastline_complexity_largest": float(
                exposed.get("coastline_complexity_largest_component", 0.0)
            ),
            "continental_ribbon_fraction_gt_0_5": float(
                continental.get("ribbon_area_fraction_gt_0_5", 0.0)
            ),
            "land_width_p50_steps": float(exposed.get("width_p50_steps", 0.0)),
            "land_width_p90_steps": float(exposed.get("width_p90_steps", 0.0)),
            "land_narrow_necks_per_1000": float(
                exposed.get("narrow_neck_cells_per_1000_mask_cells", 0.0)
            ),
            "continental_width_p50_steps": float(continental.get("width_p50_steps", 0.0)),
            "continental_width_p90_steps": float(continental.get("width_p90_steps", 0.0)),
        }
    else:
        values = {
            "land_fraction": float(
                entry_land_fraction
                if entry_land_fraction is not None
                else morphology_metrics.get("land_fraction", 0.0)
            ),
            "largest_land_component_fraction": float(
                morphology_metrics.get("largest_land_component_fraction", 0.0)
            ),
            "land_component_count": float(
                morphology_metrics.get("land_component_count", 0.0)
            ),
            "land_ribbon_fraction_gt_0_5": float(
                morphology_metrics.get("land_ribbon_fraction_gt_0_5", 0.0)
            ),
            "land_coastline_complexity_largest": float(
                morphology_metrics.get("land_coastline_complexity_largest", 0.0)
            ),
            "continental_ribbon_fraction_gt_0_5": float(
                morphology_metrics.get("continental_ribbon_fraction_gt_0_5", 0.0)
            ),
            "land_width_p50_steps": float(morphology_metrics.get("land_width_p50_steps", 0.0)),
            "land_width_p90_steps": float(morphology_metrics.get("land_width_p90_steps", 0.0)),
            "land_narrow_necks_per_1000": float(
                morphology_metrics.get("land_narrow_necks_per_1000", 0.0)
            ),
            "continental_width_p50_steps": float(
                morphology_metrics.get("continental_width_p50_steps", 0.0)
            ),
            "continental_width_p90_steps": float(
                morphology_metrics.get("continental_width_p90_steps", 0.0)
            ),
        }

    metrics: dict[str, dict[str, Any]] = {}
    out_of_envelope: list[str] = []
    for key in PLANFORM_KEYS:
        spec = REFERENCE_ENVELOPES[key]
        lo, hi = spec["range"]
        value = values[key]
        in_envelope = bool(float(lo) <= value <= float(hi))
        if not in_envelope:
            out_of_envelope.append(key)
        metrics[key] = {
            "value": value,
            "reference_min": float(lo),
            "reference_max": float(hi),
            "in_envelope": in_envelope,
            "basis": spec["basis"],
            "severity_ratio": _severity_ratio(value, float(lo), float(hi)),
        }

    flags = {
        "land_fraction_outside_initial_reference": (
            "land_fraction" in out_of_envelope
        ),
        "largest_landmass_outside_initial_reference": (
            "largest_land_component_fraction" in out_of_envelope
        ),
        "component_count_outside_initial_reference": (
            "land_component_count" in out_of_envelope
        ),
        "ribbon_outside_initial_reference": (
            "land_ribbon_fraction_gt_0_5" in out_of_envelope
        ),
        "coastline_outside_initial_reference": (
            "land_coastline_complexity_largest" in out_of_envelope
        ),
    }
    dominant_gaps = sorted(
        (
            {
                "metric": key,
                "value": metrics[key]["value"],
                "reference_min": metrics[key]["reference_min"],
                "reference_max": metrics[key]["reference_max"],
                "severity_ratio": metrics[key]["severity_ratio"],
            }
            for key in out_of_envelope
        ),
        key=lambda item: float(item["severity_ratio"]),
        reverse=True,
    )
    within_planform = not out_of_envelope
    return {
        "schema": SCHEMA,
        "status": (
            "within_initial_planform_envelope"
            if within_planform
            else "needs_planform_calibration"
        ),
        "reference_sources": REFERENCE_SOURCES,
        "metrics": metrics,
        "supporting_metrics": values,
        "out_of_envelope": out_of_envelope,
        "dominant_gaps": dominant_gaps,
        "diagnostic_flags": flags,
        "resolution_ladder": RESOLUTION_LADDER,
        "acceptance": {
            "screening_only": True,
            "source_is_direct_etopo_or_plate_model": False,
            "within_initial_planform_envelope": within_planform,
            "needs_planform_repair": bool(out_of_envelope),
            "resolution_ladder_declared": _resolution_ladder_declared(),
        },
        "notes": [
            "Use this as a tuning and regression screen before direct ETOPO/GPlates sampling.",
            "The generated Earthlike target should improve planform metrics without hiding missing tectonic parent objects.",
        ],
    }


def _severity_ratio(value: float, reference_min: float, reference_max: float) -> float:
    if value < reference_min:
        return float(reference_min / max(value, 1.0e-12))
    if value > reference_max:
        return float(value / max(reference_max, 1.0e-12))
    return 1.0


def _resolution_ladder_declared() -> bool:
    cells = {str(item.get("cells", "")) for item in RESOLUTION_LADDER}
    return {"8000", "24000", "72000"} <= cells
