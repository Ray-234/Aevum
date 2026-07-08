"""Boundary-process geometry references for Stage-B tectonics diagnostics."""
from __future__ import annotations

from collections import Counter, deque
import hashlib
import json
from typing import Any

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.physiographic_reference import REFERENCE_SOURCES


SCHEMA = "aevum.boundary_process_geometry.v1"

SOURCE_IDS = (
    "PB2002_PLATE_BOUNDARIES",
    "GPLATES",
    "PYGPLATES_EARTHBYTE",
    "EARTHBYTE_RECONSTRUCTIONS",
    "GEM_GLOBAL_ACTIVE_FAULTS",
)

PROCESS_TYPES = (
    "ridge",
    "transform",
    "subduction_trench",
    "collision_suture",
    "diffuse_deformation",
    "passive_margin",
    "continental_rift",
)

CORE_PROCESS_TYPES = (
    "ridge",
    "transform",
    "subduction_trench",
    "collision_suture",
    "diffuse_deformation",
)

EXPECTED_CURRENT_RESIDUAL_TYPES: tuple[str, ...] = ()

REFERENCE_LENGTH_FRACTION_ENVELOPES = {
    "ridge": {"min": 0.08, "max": 0.35},
    "transform": {"min": 0.05, "max": 0.25},
    "subduction_trench": {"min": 0.08, "max": 0.35},
    "collision_suture": {"min": 0.02, "max": 0.18},
    "diffuse_deformation": {"min": 0.04, "max": 0.30},
    "passive_margin": {"min": 0.05, "max": 0.35},
    "continental_rift": {"min": 0.02, "max": 0.20},
}


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def boundary_process_geometry_summary(world=None, *, cells: int = 1200) -> dict[str, Any]:
    """Return reference, synthetic, and optional generated-world boundary metrics."""
    reference = boundary_process_geometry_reference_summary(cells=cells)
    current = None if world is None else generated_world_boundary_geometry_summary(world)
    acceptance = {
        "reference_geometry_ready": reference["status"] == "boundary_process_geometry_reference_ready",
        "generated_world_comparison_available": current is not None,
    }
    if current is not None:
        acceptance.update({
            "current_boundary_network_present": current["acceptance"]["boundary_network_present"],
            "current_key_types_present": current["acceptance"]["key_boundary_types_present"],
            "current_expected_residuals_recorded": current["acceptance"]["expected_residuals_recorded"],
            "current_unexpected_missing_types_absent": current["acceptance"]["unexpected_missing_types_absent"],
        })
    summary = {
        "schema": SCHEMA,
        "status": (
            "boundary_process_geometry_ready"
            if all(acceptance.values())
            else "boundary_process_geometry_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "next_gates": (
            "P82.wilson_cycle_lifecycle_reference",
            "P83.crust_sediment_province_coupling",
        ),
    }
    summary["summary_digest"] = _digest(summary)
    return summary


def boundary_process_geometry_reference_summary(*, cells: int = 1200) -> dict[str, Any]:
    """Small derived reference and synthetic network fixture for P81."""
    known_source_ids = {str(source["id"]) for source in REFERENCE_SOURCES}
    missing_source_ids = tuple(sorted(set(SOURCE_IDS) - known_source_ids))
    grid, raw_network = _synthetic_boundary_network(cells=cells)
    synthetic = analyze_boundary_network(
        grid,
        raw_network,
        source="synthetic_reference_boundary_network",
    )
    missing_synthetic_types = tuple(
        sorted(set(PROCESS_TYPES) - set(synthetic["present_process_types"]))
    )
    acceptance = {
        "source_ids_known": not missing_source_ids,
        "reference_envelopes_complete": set(REFERENCE_LENGTH_FRACTION_ENVELOPES) == set(PROCESS_TYPES),
        "small_derived_fixture": True,
        "raw_vectors_not_stored": True,
        "direct_vector_extraction_marked_pending": True,
        "synthetic_required_types_covered": not missing_synthetic_types,
        "synthetic_type_diversity": synthetic["process_type_count"] >= 7,
        "synthetic_length_fractions_in_reference_envelope": all(
            check["in_envelope"] for check in synthetic["length_fraction_checks"].values()
        ),
        "synthetic_spherical_continuity": synthetic["antimeridian_ridge_component_count"] == 1,
        "synthetic_ridge_transform_adjacency": (
            synthetic["adjacency"]["transform_near_ridge_fraction"] >= 0.65
        ),
        "synthetic_transform_offsets": synthetic["transform_offset_count"] >= 2,
        "synthetic_trench_active_margin_adjacency": (
            synthetic["adjacency"]["trench_near_active_margin_fraction"] >= 0.75
        ),
        "synthetic_collision_diffuse_adjacency": (
            synthetic["adjacency"]["collision_near_diffuse_fraction"] >= 0.60
        ),
        "synthetic_geometry_not_single_straight_line": (
            synthetic["mean_component_sinuosity"] >= 1.05
            and synthetic["trench_longitude_std_deg"] >= 2.0
        ),
    }
    summary = {
        "schema": "aevum.boundary_process_geometry_reference.v1",
        "source_ids": SOURCE_IDS,
        "missing_source_ids": missing_source_ids,
        "process_types": PROCESS_TYPES,
        "core_process_types": CORE_PROCESS_TYPES,
        "reference_length_fraction_envelopes": REFERENCE_LENGTH_FRACTION_ENVELOPES,
        "synthetic_network": synthetic,
        "missing_synthetic_types": missing_synthetic_types,
        "extraction_policy": {
            "derived_from": "small deterministic spherical process-network fixture",
            "raw_pb2002_or_gplates_vectors_stored": False,
            "direct_vector_extraction_pending": True,
            "exact_earth_geometry_required": False,
            "process_geometry_required": True,
        },
        "acceptance": acceptance,
        "next_gates": (
            "P82.wilson_cycle_lifecycle_reference",
            "P83.crust_sediment_province_coupling",
        ),
    }
    summary["status"] = (
        "boundary_process_geometry_reference_ready"
        if all(acceptance.values())
        else "boundary_process_geometry_reference_incomplete"
    )
    summary["fixture_digest"] = _digest(summary)
    return summary


def generated_world_boundary_geometry_summary(world) -> dict[str, Any]:
    """Analyze a generated world's current boundary network with P81 schema."""
    raw_boundaries = world.networks.get("tectonics.boundaries", {})
    network = {
        "ridge": raw_boundaries.get("ridge", ()),
        "transform": raw_boundaries.get("transform", ()),
        "subduction_trench": _merge_cells(
            raw_boundaries.get("trench", ()),
            raw_boundaries.get("subduction", ()),
        ),
        "collision_suture": _merge_cells(
            raw_boundaries.get("collision", ()),
            raw_boundaries.get("suture", ()),
        ),
        "diffuse_deformation": _merge_cells(
            raw_boundaries.get("convergent", ()),
            raw_boundaries.get("divergent", ()),
        ),
        "passive_margin": raw_boundaries.get("passive_margin", ()),
        "continental_rift": raw_boundaries.get("divergent", ()),
        "active_margin": raw_boundaries.get("active_margin", ()),
    }
    metrics = analyze_boundary_network(
        world.grid,
        network,
        source="generated_world_tectonics_boundaries",
    )
    present = set(metrics["present_process_types"])
    missing_types = tuple(sorted(set(PROCESS_TYPES) - present))
    unexpected_missing = tuple(sorted(set(missing_types) - set(EXPECTED_CURRENT_RESIDUAL_TYPES)))
    acceptance = {
        "boundary_network_present": metrics["any_boundary_cell_count"] > 0,
        "key_boundary_types_present": not unexpected_missing,
        "expected_residuals_recorded": (
            set(missing_types).issubset(set(EXPECTED_CURRENT_RESIDUAL_TYPES))
        ),
        "unexpected_missing_types_absent": not unexpected_missing,
        "ridge_and_trench_present": (
            "ridge" in present and "subduction_trench" in present
        ),
        "rift_or_passive_margin_present": (
            "continental_rift" in present or "passive_margin" in present
        ),
        "collision_or_suture_present": "collision_suture" in present,
    }
    summary = {
        "schema": "aevum.generated_boundary_process_geometry.v1",
        "context": {
            "spec_name": str(world.spec.name),
            "seed": int(world.spec.seed),
            "cells": int(world.grid.n),
            "time_myr": float(world.time_myr),
        },
        "metrics": metrics,
        "missing_process_types": missing_types,
        "unexpected_missing_process_types": unexpected_missing,
        "expected_current_residual_types": EXPECTED_CURRENT_RESIDUAL_TYPES,
        "acceptance": acceptance,
        "limitations": {
            "current_transform_boundary_missing": "transform" in missing_types,
            "current_boundary_network_is_cell_mask_not_vector_topology": True,
            "direct_gplates_vector_comparison_pending": True,
        },
    }
    summary["status"] = (
        "generated_boundary_process_geometry_ready"
        if all(acceptance.values())
        else "generated_boundary_process_geometry_incomplete"
    )
    summary["summary_digest"] = _digest(summary)
    return summary


def analyze_boundary_network(
    grid: SphereGrid,
    raw_network: dict[str, Any],
    *,
    source: str,
) -> dict[str, Any]:
    masks = _normalise_network(grid, raw_network)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    length_weight = np.sqrt(np.maximum(area, 1.0))
    type_lengths = {
        name: float(length_weight[mask].sum())
        for name, mask in masks.items()
        if name in PROCESS_TYPES
    }
    total_length = max(sum(type_lengths.values()), 1.0e-12)
    length_fractions = {
        name: float(type_lengths.get(name, 0.0) / total_length)
        for name in PROCESS_TYPES
    }
    length_fraction_checks = {
        name: {
            "value": float(length_fractions[name]),
            "min": float(REFERENCE_LENGTH_FRACTION_ENVELOPES[name]["min"]),
            "max": float(REFERENCE_LENGTH_FRACTION_ENVELOPES[name]["max"]),
            "in_envelope": (
                REFERENCE_LENGTH_FRACTION_ENVELOPES[name]["min"]
                <= length_fractions[name]
                <= REFERENCE_LENGTH_FRACTION_ENVELOPES[name]["max"]
            ),
        }
        for name in PROCESS_TYPES
    }
    type_metrics = {
        name: _mask_geometry_metrics(grid, mask)
        for name, mask in masks.items()
        if name in PROCESS_TYPES
    }
    all_mask = np.zeros(grid.n, dtype=bool)
    for name in PROCESS_TYPES:
        all_mask |= masks[name]
    present = tuple(name for name in PROCESS_TYPES if masks[name].any())
    component_sinuosities = [
        value
        for stats in type_metrics.values()
        for value in stats["component_sinuosities"]
    ]
    trench_lons = grid.lon[masks["subduction_trench"]]
    summary = {
        "schema": "aevum.boundary_network_geometry_metrics.v1",
        "source": str(source),
        "cells": int(grid.n),
        "process_types": PROCESS_TYPES,
        "present_process_types": present,
        "process_type_count": len(present),
        "any_boundary_cell_count": int(all_mask.sum()),
        "boundary_cell_fraction": float(all_mask.sum() / max(grid.n, 1)),
        "type_cell_counts": {
            name: int(masks[name].sum()) for name in PROCESS_TYPES
        },
        "type_length_proxy_m": type_lengths,
        "type_length_fractions": length_fractions,
        "length_fraction_checks": length_fraction_checks,
        "type_geometry": type_metrics,
        "mean_component_sinuosity": (
            float(np.mean(component_sinuosities)) if component_sinuosities else 0.0
        ),
        "max_component_sinuosity": (
            float(np.max(component_sinuosities)) if component_sinuosities else 0.0
        ),
        "trench_longitude_std_deg": (
            float(np.std(trench_lons)) if trench_lons.size else 0.0
        ),
        "transform_offset_count": _component_count(grid, masks["transform"]),
        "antimeridian_ridge_component_count": _component_count(
            grid,
            masks.get("antimeridian_ridge", np.zeros(grid.n, dtype=bool)),
        ),
        "adjacency": {
            "transform_near_ridge_fraction": _near_fraction(
                grid, masks["transform"], masks["ridge"], passes=1),
            "ridge_near_transform_fraction": _near_fraction(
                grid, masks["ridge"], masks["transform"], passes=1),
            "trench_near_active_margin_fraction": _near_fraction(
                grid, masks["subduction_trench"], masks["active_margin"], passes=2),
            "collision_near_diffuse_fraction": _near_fraction(
                grid, masks["collision_suture"], masks["diffuse_deformation"], passes=2),
            "rift_near_passive_margin_fraction": _near_fraction(
                grid, masks["continental_rift"], masks["passive_margin"], passes=3),
        },
    }
    return summary


def _synthetic_boundary_network(*, cells: int) -> tuple[SphereGrid, dict[str, np.ndarray]]:
    grid = SphereGrid.fibonacci(cells, CONSTANTS.EARTH_RADIUS)
    lat = grid.lat
    lon = grid.lon
    ridge = (
        ((np.abs(lon + 40.0) < 4.0) & (lat > -66.0) & (lat <= -24.0))
        | ((np.abs(lon + 15.0) < 4.0) & (lat > -24.0) & (lat <= 15.0))
        | ((np.abs(lon - 10.0) < 4.0) & (lat > 15.0) & (lat < 66.0))
    )
    transform = (
        ((np.abs(lat + 24.0) < 4.0) & (lon > -43.0) & (lon < -12.0))
        | ((np.abs(lat - 15.0) < 4.0) & (lon > -18.0) & (lon < 13.0))
    )
    trench_center = 121.0 + 11.0 * np.sin(np.radians(2.3 * lat))
    trench = (np.abs(lon - trench_center) < 4.0) & (lat > -58.0) & (lat < 58.0)
    active_margin = (
        (np.abs(lon - (trench_center + 7.5)) < 5.0)
        & (lat > -58.0)
        & (lat < 58.0)
    )
    collision_axis = 7.0 * np.sin(np.radians(2.0 * (lon - 70.0)))
    collision = (
        (np.abs(lat - collision_axis) < 4.0)
        & (lon > 34.0)
        & (lon < 103.0)
    )
    diffuse = (
        (np.abs(lat - collision_axis) < 9.0)
        & (lon > 29.0)
        & (lon < 108.0)
        & ~collision
    )
    passive_margin = (np.abs(lon + 128.0) < 4.0) & (lat > -55.0) & (lat < 55.0)
    continental_rift = (
        (np.abs(lon - 70.0) < 4.0)
        & (lat > -53.0)
        & (lat < -8.0)
    )
    antimeridian_ridge = (
        (np.abs(np.abs(lon) - 178.0) < 5.5)
        & (lat > -23.0)
        & (lat < 23.0)
    )
    ridge = ridge | antimeridian_ridge
    return grid, {
        "ridge": np.where(ridge)[0],
        "transform": np.where(transform)[0],
        "subduction_trench": np.where(trench)[0],
        "active_margin": np.where(active_margin | trench)[0],
        "collision_suture": np.where(collision)[0],
        "diffuse_deformation": np.where(diffuse)[0],
        "passive_margin": np.where(passive_margin)[0],
        "continental_rift": np.where(continental_rift)[0],
        "antimeridian_ridge": np.where(antimeridian_ridge)[0],
    }


def _normalise_network(grid: SphereGrid, raw_network: dict[str, Any]) -> dict[str, np.ndarray]:
    masks = {name: np.zeros(grid.n, dtype=bool) for name in PROCESS_TYPES}
    masks["active_margin"] = np.zeros(grid.n, dtype=bool)
    masks["antimeridian_ridge"] = np.zeros(grid.n, dtype=bool)
    aliases = {
        "ridge": ("ridge",),
        "transform": ("transform",),
        "subduction_trench": ("subduction_trench", "trench", "subduction"),
        "collision_suture": ("collision_suture", "collision", "suture"),
        "diffuse_deformation": ("diffuse_deformation", "convergent", "divergent"),
        "passive_margin": ("passive_margin",),
        "continental_rift": ("continental_rift", "rift", "divergent"),
        "active_margin": ("active_margin",),
        "antimeridian_ridge": ("antimeridian_ridge",),
    }
    for canonical, names in aliases.items():
        cells: list[int] = []
        for name in names:
            cells.extend(int(c) for c in np.asarray(raw_network.get(name, ()), dtype=int).ravel())
        if not cells:
            continue
        arr = np.asarray(cells, dtype=int)
        arr = arr[(arr >= 0) & (arr < grid.n)]
        masks[canonical][np.unique(arr)] = True
    return masks


def _mask_geometry_metrics(grid: SphereGrid, mask: np.ndarray) -> dict[str, Any]:
    components = _components(grid, mask)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    length_weight = np.sqrt(np.maximum(area, 1.0))
    lengths = [float(length_weight[comp].sum()) for comp in components]
    sinuosities = [_component_sinuosity(grid, comp) for comp in components]
    largest = max(lengths, default=0.0)
    total = max(sum(lengths), 1.0e-12)
    return {
        "cell_count": int(mask.sum()),
        "component_count": len(components),
        "largest_component_length_fraction": float(largest / total),
        "endpoint_count": _endpoint_count(grid, mask),
        "component_sinuosities": tuple(round(float(v), 4) for v in sinuosities),
        "mean_sinuosity": float(np.mean(sinuosities)) if sinuosities else 0.0,
        "max_sinuosity": float(np.max(sinuosities)) if sinuosities else 0.0,
    }


def _components(grid: SphereGrid, mask: np.ndarray) -> list[np.ndarray]:
    mask = np.asarray(mask, dtype=bool)
    seen = np.zeros(grid.n, dtype=bool)
    out: list[np.ndarray] = []
    for start in np.where(mask)[0]:
        start = int(start)
        if seen[start]:
            continue
        cells: list[int] = []
        queue: deque[int] = deque([start])
        seen[start] = True
        while queue:
            c = queue.popleft()
            cells.append(c)
            for nb in grid.neighbors[c]:
                nb = int(nb)
                if mask[nb] and not seen[nb]:
                    seen[nb] = True
                    queue.append(nb)
        out.append(np.asarray(cells, dtype=int))
    return out


def _component_count(grid: SphereGrid, mask: np.ndarray) -> int:
    return len(_components(grid, np.asarray(mask, dtype=bool)))


def _component_sinuosity(grid: SphereGrid, cells: np.ndarray) -> float:
    cells = np.asarray(cells, dtype=int)
    if cells.size <= 1:
        return 1.0
    area = np.asarray(grid.cell_area, dtype=np.float64)
    length = float(np.sqrt(np.maximum(area[cells], 1.0)).sum())
    xyz = grid.xyz[cells]
    dots = np.clip(xyz @ xyz.T, -1.0, 1.0)
    diameter = float(np.max(np.arccos(dots)) * grid.radius_m)
    cell_scale = float(np.sqrt(np.mean(area[cells])))
    return float(length / max(diameter, cell_scale, 1.0))


def _endpoint_count(grid: SphereGrid, mask: np.ndarray) -> int:
    mask = np.asarray(mask, dtype=bool)
    count = 0
    for c in np.where(mask)[0]:
        n_same = int(np.count_nonzero(mask[grid.neighbors[int(c)]]))
        if n_same <= 1:
            count += 1
    return count


def _dilate(grid: SphereGrid, mask: np.ndarray, *, passes: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(0, int(passes))):
        nxt = out.copy()
        for c in np.where(out)[0]:
            nxt[grid.neighbors[int(c)]] = True
        out = nxt
    return out


def _near_fraction(
    grid: SphereGrid,
    source: np.ndarray,
    target: np.ndarray,
    *,
    passes: int,
) -> float:
    source = np.asarray(source, dtype=bool)
    target = np.asarray(target, dtype=bool)
    if not source.any():
        return 0.0
    near = _dilate(grid, target, passes=passes)
    return float(np.mean(near[source]))


def _merge_cells(*arrays: Any) -> np.ndarray:
    cells: list[int] = []
    for arr in arrays:
        cells.extend(int(c) for c in np.asarray(arr, dtype=int).ravel())
    if not cells:
        return np.array([], dtype=int)
    return np.unique(np.asarray(cells, dtype=int))
