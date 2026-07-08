"""Drainage divide and province-alignment reference diagnostics for P85."""
from __future__ import annotations

import hashlib
import json
from collections import Counter, deque
from typing import Any

import numpy as np


SCHEMA = "aevum.drainage_divide_province_alignment.v1"

SOURCE_IDS = ("HYDROSHEDS_HYDROBASINS_HYDRORIVERS",)

DIVIDE = "D"

EXPECTED_CURRENT_RESIDUAL_ITEMS = (
    "terrain.drainage_basins",
    "terrain.drainage_divides",
    "terrain.flow_direction",
    "terrain.flow_accumulation",
)

REQUIRED_PROVINCE_CLASSES = {
    "shield",
    "platform",
    "intracratonic_basin",
    "active_orogen",
    "foreland_basin",
    "rift_shoulder",
    "rift_axis",
    "rift_basin",
    "passive_margin_lowland",
    "continental_shelf",
}

SINK_CLASSES = {
    "west_interior": {"intracratonic_basin"},
    "east_margin": {"passive_margin_lowland", "continental_shelf"},
    "south_rift": {"rift_basin", "rift_axis"},
}

DIVIDE_PROVINCE_CLASSES = {
    "active_orogen",
    "rift_shoulder",
    "rift_axis",
    "old_suture",
}


def _rows(*rows: str) -> tuple[tuple[str, ...], ...]:
    return tuple(tuple(row.split()) for row in rows)


REFERENCE_FIXTURE: dict[str, Any] = {
    "name": "p85_multi_province_drainage_reference",
    "source_ids": SOURCE_IDS,
    "province_grid": _rows(
        "SH SH PL AO FB FB PM PM PM CS",
        "SH SH PL AO FB FB PM PM PM CS",
        "SH PL PL AO FB FB PM PM PM CS",
        "IB IB PL AO FB FB PM PM CS CS",
        "IB IB PL AO FB FB PM PM CS CS",
        "PL PL RS RS RB RB RB PM CS CS",
        "PL PL RS RA RB RB RB PM CS CS",
        "PL PL RS RA RB RB RB PM CS CS",
    ),
    "basin_grid": _rows(
        "W W W D E E E E E E",
        "W W W D E E E E E E",
        "W W W D E E E E E E",
        "W W W D E E E E E E",
        "W W W D E E E E E E",
        "R R R R R R R E E E",
        "R R R R R R R R E E",
        "R R R R R R R R E E",
    ),
    "province_classes": {
        "SH": "shield",
        "PL": "platform",
        "IB": "intracratonic_basin",
        "AO": "active_orogen",
        "FB": "foreland_basin",
        "RS": "rift_shoulder",
        "RA": "rift_axis",
        "RB": "rift_basin",
        "PM": "passive_margin_lowland",
        "CS": "continental_shelf",
    },
    "basin_ids": {
        "W": "west_interior",
        "E": "east_margin",
        "R": "south_rift",
        "D": "divide",
    },
    "elevation_m": {
        "SH": 680.0,
        "PL": 320.0,
        "IB": 90.0,
        "AO": 2850.0,
        "FB": 180.0,
        "RS": 1150.0,
        "RA": 120.0,
        "RB": -40.0,
        "PM": 35.0,
        "CS": -100.0,
    },
    "flow_paths": (
        {
            "id": "west-shield-to-interior-basin",
            "basin_id": "west_interior",
            "cells": ((0, 0), (1, 1), (2, 2), (3, 1)),
            "expected_sink_class": "intracratonic_basin",
        },
        {
            "id": "west-platform-to-interior-basin",
            "basin_id": "west_interior",
            "cells": ((2, 1), (3, 2), (4, 1)),
            "expected_sink_class": "intracratonic_basin",
        },
        {
            "id": "orogen-foreland-to-passive-margin",
            "basin_id": "east_margin",
            "cells": ((0, 4), (1, 5), (2, 6), (3, 7), (4, 8), (4, 9)),
            "expected_sink_class": "continental_shelf",
        },
        {
            "id": "foreland-to-coastal-plain",
            "basin_id": "east_margin",
            "cells": ((3, 4), (3, 5), (4, 6), (5, 7), (5, 8)),
            "expected_sink_class": "continental_shelf",
        },
        {
            "id": "rift-shoulder-to-rift-basin",
            "basin_id": "south_rift",
            "cells": ((5, 2), (5, 3), (5, 4), (6, 5), (7, 6)),
            "expected_sink_class": "rift_basin",
        },
        {
            "id": "southern-rift-shoulder-to-rift-axis",
            "basin_id": "south_rift",
            "cells": ((7, 2), (7, 3), (7, 4)),
            "expected_sink_class": "rift_basin",
        },
    ),
}


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _grid_shape(grid: tuple[tuple[str, ...], ...]) -> tuple[int, int]:
    return len(grid), len(grid[0]) if grid else 0


def _neighbors(r: int, c: int, rows: int, cols: int) -> tuple[tuple[int, int], ...]:
    out = []
    for nr, nc in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
        if 0 <= nr < rows and 0 <= nc < cols:
            out.append((nr, nc))
    return tuple(out)


def _province_class(fixture: dict[str, Any], r: int, c: int) -> str:
    province_code = fixture["province_grid"][r][c]
    return str(fixture["province_classes"][province_code])


def _elevation(fixture: dict[str, Any], r: int, c: int) -> float:
    province_code = fixture["province_grid"][r][c]
    return float(fixture["elevation_m"][province_code])


def _component_counts(
    basin_grid: tuple[tuple[str, ...], ...],
) -> dict[str, int]:
    rows, cols = _grid_shape(basin_grid)
    seen: set[tuple[int, int]] = set()
    components: Counter[str] = Counter()
    for r in range(rows):
        for c in range(cols):
            basin_code = basin_grid[r][c]
            if basin_code == DIVIDE or (r, c) in seen:
                continue
            components[basin_code] += 1
            queue: deque[tuple[int, int]] = deque([(r, c)])
            seen.add((r, c))
            while queue:
                cr, cc = queue.popleft()
                for nr, nc in _neighbors(cr, cc, rows, cols):
                    if (nr, nc) in seen:
                        continue
                    if basin_grid[nr][nc] == basin_code:
                        seen.add((nr, nc))
                        queue.append((nr, nc))
    return {
        basin_code: int(count)
        for basin_code, count in sorted(components.items())
    }


def _divide_alignment_summary(fixture: dict[str, Any]) -> dict[str, Any]:
    province_grid = fixture["province_grid"]
    basin_grid = fixture["basin_grid"]
    rows, cols = _grid_shape(basin_grid)
    divide_cells = [
        (r, c)
        for r in range(rows)
        for c in range(cols)
        if basin_grid[r][c] == DIVIDE
    ]
    aligned = []
    boundary_aligned = []
    highland_aligned = []
    for r, c in divide_cells:
        province_class = _province_class(fixture, r, c)
        adjacent_basin_codes = {
            basin_grid[nr][nc]
            for nr, nc in _neighbors(r, c, rows, cols)
            if basin_grid[nr][nc] != DIVIDE
        }
        adjacent_province_classes = {
            _province_class(fixture, nr, nc)
            for nr, nc in _neighbors(r, c, rows, cols)
            if province_grid[nr][nc] != province_grid[r][c]
        }
        is_highland = province_class in DIVIDE_PROVINCE_CLASSES
        is_boundary = len(adjacent_basin_codes) >= 2 or bool(adjacent_province_classes)
        highland_aligned.append(is_highland)
        boundary_aligned.append(is_boundary)
        aligned.append(is_highland or is_boundary)
    divide_count = len(divide_cells)
    return {
        "divide_cell_count": divide_count,
        "aligned_divide_cell_count": int(sum(aligned)),
        "highland_aligned_divide_cell_count": int(sum(highland_aligned)),
        "boundary_aligned_divide_cell_count": int(sum(boundary_aligned)),
        "divide_alignment_fraction": (
            float(sum(aligned) / divide_count) if divide_count else 0.0
        ),
        "highland_alignment_fraction": (
            float(sum(highland_aligned) / divide_count) if divide_count else 0.0
        ),
        "boundary_alignment_fraction": (
            float(sum(boundary_aligned) / divide_count) if divide_count else 0.0
        ),
    }


def _flow_path_summary(fixture: dict[str, Any]) -> dict[str, Any]:
    basin_grid = fixture["basin_grid"]
    basin_ids = fixture["basin_ids"]
    paths = fixture["flow_paths"]
    path_summaries = []
    sink_failures = []
    basin_crossings = []
    divide_crossings = []
    uphill_steps = 0
    total_steps = 0
    for path in paths:
        cells = tuple(tuple(cell) for cell in path["cells"])
        basin_id = str(path["basin_id"])
        sink_r, sink_c = cells[-1]
        sink_class = _province_class(fixture, sink_r, sink_c)
        basin_codes = tuple(basin_grid[r][c] for r, c in cells)
        resolved_basin_ids = tuple(str(basin_ids[code]) for code in basin_codes)
        elevations = tuple(_elevation(fixture, r, c) for r, c in cells)
        path_uphill_steps = sum(
            1
            for left, right in zip(elevations, elevations[1:])
            if right > left + 1.0e-9
        )
        total_steps += max(len(cells) - 1, 0)
        uphill_steps += path_uphill_steps
        crossed_divide = DIVIDE in basin_codes
        crossed_other_basin = any(
            resolved not in {basin_id, "divide"} for resolved in resolved_basin_ids)
        sink_ok = sink_class in SINK_CLASSES[basin_id] and (
            sink_class == str(path["expected_sink_class"])
            or sink_class in SINK_CLASSES[basin_id]
        )
        if crossed_divide:
            divide_crossings.append(str(path["id"]))
        if crossed_other_basin:
            basin_crossings.append(str(path["id"]))
        if not sink_ok:
            sink_failures.append(str(path["id"]))
        path_summaries.append({
            "id": str(path["id"]),
            "basin_id": basin_id,
            "sink_class": sink_class,
            "expected_sink_class": str(path["expected_sink_class"]),
            "step_count": int(max(len(cells) - 1, 0)),
            "uphill_step_count": int(path_uphill_steps),
            "crossed_divide": bool(crossed_divide),
            "crossed_other_basin": bool(crossed_other_basin),
            "sink_ok": bool(sink_ok),
            "elevation_drop_m": float(elevations[0] - elevations[-1]),
        })
    return {
        "path_count": int(len(paths)),
        "path_summaries": tuple(path_summaries),
        "sink_failure_count": int(len(sink_failures)),
        "sink_failures": tuple(sink_failures),
        "basin_crossing_count": int(len(basin_crossings)),
        "basin_crossings": tuple(basin_crossings),
        "divide_crossing_count": int(len(divide_crossings)),
        "divide_crossings": tuple(divide_crossings),
        "uphill_step_count": int(uphill_steps),
        "total_step_count": int(total_steps),
        "flow_to_sink_consistency_fraction": (
            float((len(paths) - len(sink_failures) - len(basin_crossings)
                   - len(divide_crossings)) / len(paths))
            if paths else 0.0
        ),
        "downhill_step_fraction": (
            float((total_steps - uphill_steps) / total_steps)
            if total_steps else 1.0
        ),
    }


def reference_drainage_divide_summary() -> dict[str, Any]:
    fixture = REFERENCE_FIXTURE
    province_grid = fixture["province_grid"]
    basin_grid = fixture["basin_grid"]
    rows, cols = _grid_shape(province_grid)
    basin_rows, basin_cols = _grid_shape(basin_grid)
    classes_present = sorted({
        _province_class(fixture, r, c)
        for r in range(rows)
        for c in range(cols)
    })
    missing_required_classes = tuple(
        sorted(REQUIRED_PROVINCE_CLASSES - set(classes_present)))
    basin_codes = {
        basin_grid[r][c]
        for r in range(basin_rows)
        for c in range(basin_cols)
        if basin_grid[r][c] != DIVIDE
    }
    component_counts = _component_counts(basin_grid)
    max_component_count = max(component_counts.values()) if component_counts else 0
    divide = _divide_alignment_summary(fixture)
    flow = _flow_path_summary(fixture)
    total_cells = basin_rows * basin_cols
    divide_fraction = (
        float(divide["divide_cell_count"] / total_cells) if total_cells else 0.0
    )
    checkerboard_free = bool(
        max_component_count == 1 and 0.05 <= divide_fraction <= 0.20)
    acceptance = {
        "fixture_schema_ready": bool(fixture["name"] and SOURCE_IDS),
        "grid_shapes_match": bool((rows, cols) == (basin_rows, basin_cols)),
        "required_province_classes_present": not missing_required_classes,
        "required_basin_ids_present": {"W", "E", "R"}.issubset(basin_codes),
        "divide_boundary_alignment": divide["divide_alignment_fraction"] >= 0.90,
        "divide_highland_alignment": divide["highland_alignment_fraction"] >= 0.75,
        "flow_to_expected_sinks": flow["flow_to_sink_consistency_fraction"] >= 0.95,
        "flow_paths_do_not_cross_divides": flow["divide_crossing_count"] == 0,
        "flow_paths_stay_in_basins": flow["basin_crossing_count"] == 0,
        "flow_paths_are_downhill": flow["uphill_step_count"] == 0,
        "basins_not_checkerboarded": checkerboard_free,
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "row_count": int(rows),
        "col_count": int(cols),
        "province_class_count": int(len(classes_present)),
        "basin_count": int(len(basin_codes)),
        "divide_cell_count": int(divide["divide_cell_count"]),
        "divide_fraction": float(divide_fraction),
        "divide_alignment_fraction": float(divide["divide_alignment_fraction"]),
        "highland_alignment_fraction": float(divide["highland_alignment_fraction"]),
        "boundary_alignment_fraction": float(divide["boundary_alignment_fraction"]),
        "flow_path_count": int(flow["path_count"]),
        "flow_to_sink_consistency_fraction": float(
            flow["flow_to_sink_consistency_fraction"]),
        "downhill_step_fraction": float(flow["downhill_step_fraction"]),
        "uphill_step_count": int(flow["uphill_step_count"]),
        "divide_crossing_count": int(flow["divide_crossing_count"]),
        "basin_crossing_count": int(flow["basin_crossing_count"]),
        "sink_failure_count": int(flow["sink_failure_count"]),
        "max_basin_component_count": int(max_component_count),
        "missing_required_province_class_count": int(len(missing_required_classes)),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "drainage_divide_reference_ready"
            if all(acceptance.values())
            else "drainage_divide_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "fixture_name": fixture["name"],
        "classes_present": tuple(classes_present),
        "missing_required_province_classes": missing_required_classes,
        "basin_component_counts": component_counts,
        "divide": divide,
        "flow": flow,
        "extraction_policy": {
            "raw_hydrology_data_stored": False,
            "direct_hydrology_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest(fixture),
    }


def current_generated_drainage_audit(world: Any) -> dict[str, Any]:
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    networks = getattr(world, "networks", {})
    missing_items = []
    for name in EXPECTED_CURRENT_RESIDUAL_ITEMS:
        if name not in fields and name not in objects and name not in networks:
            missing_items.append(name)
    landforms = list(objects.get("terrain.continental_landforms", []))
    kind_counts = Counter(str(obj.get("kind", "unknown")) for obj in landforms)
    province_objects = list(objects.get("tectonics.continental_provinces", []))
    province_kind_counts = Counter(
        str(obj.get("province_class", "unknown")) for obj in province_objects)
    context_groups = {
        "divide_parent_context": (
            kind_counts.get("old_subdued_orogen", 0)
            + kind_counts.get("orogenic_plateau", 0)
            + kind_counts.get("active_orogen", 0)
        ),
        "basin_sink_context": (
            kind_counts.get("interior_basin", 0)
            + kind_counts.get("foreland_basin", 0)
            + kind_counts.get("rift_basin", 0)
        ),
        "margin_sink_context": kind_counts.get("passive_margin_lowland", 0),
    }
    production_context_groups = {
        "divide_parent_context": (
            province_kind_counts.get("old_orogen", 0)
            + province_kind_counts.get("active_orogen", 0)
            + province_kind_counts.get("rift_system", 0)
        ),
        "basin_sink_context": (
            province_kind_counts.get("intracratonic_basin", 0)
            + province_kind_counts.get("foreland_basin", 0)
            + province_kind_counts.get("rift_system", 0)
        ),
        "margin_sink_context": province_kind_counts.get("passive_margin_lowland", 0),
    }
    expected_residuals_recorded = bool(
        set(missing_items).issubset(set(EXPECTED_CURRENT_RESIDUAL_ITEMS)))
    n = int(world.grid.n)
    area = np.asarray(getattr(world.grid, "cell_area", np.ones(n)), dtype=np.float64)
    elevation = (
        np.asarray(fields["terrain.elevation_m"], dtype=np.float64)
        if "terrain.elevation_m" in fields else np.zeros(n, dtype=np.float64)
    )
    flow_surface = (
        np.asarray(fields["terrain.drainage_surface_m"], dtype=np.float64)
        if "terrain.drainage_surface_m" in fields else elevation
    )
    sea_level = float(getattr(world, "sea_level", 0.0))
    land = elevation >= sea_level
    basin = (
        np.asarray(fields.get("terrain.drainage_basins", np.full(n, -1.0)),
                   dtype=np.float64).astype(int)
    )
    divides = (
        np.asarray(fields.get("terrain.drainage_divides", np.zeros(n)),
                   dtype=np.float64) > 0.5
    )
    receiver = (
        np.asarray(fields.get("terrain.flow_direction", np.arange(n)),
                   dtype=np.float64).astype(int)
    )
    receiver[(receiver < 0) | (receiver >= n)] = np.arange(n)[
        (receiver < 0) | (receiver >= n)]
    province = (
        np.asarray(fields.get("terrain.continental_province_code", np.zeros(n)),
                   dtype=np.float64).astype(int)
    )
    valid = land & (basin >= 0)
    basin_ids = np.unique(basin[valid])

    component_failures = 0
    major_component_count = 0
    for bid in basin_ids:
        mask = valid & (basin == int(bid))
        basin_area = float(area[mask].sum())
        if basin_area <= 0.005 * max(float(area[valid].sum()), 1.0):
            continue
        major_component_count += 1
        seen = np.zeros(n, dtype=bool)
        components = 0
        for start in np.where(mask)[0]:
            if seen[int(start)]:
                continue
            components += 1
            stack = [int(start)]
            seen[int(start)] = True
            while stack:
                c = stack.pop()
                for nb in world.grid.neighbors[c]:
                    nb = int(nb)
                    if mask[nb] and not seen[nb]:
                        seen[nb] = True
                        stack.append(nb)
        if components > 1:
            component_failures += 1

    divide_cells = divides & valid
    boundary_aligned = np.zeros(n, dtype=bool)
    for c in np.where(divide_cells)[0]:
        c = int(c)
        nb = world.grid.neighbors[c]
        boundary_aligned[c] = bool(np.any((basin[nb] >= 0) & (basin[nb] != basin[c])))
    if valid.any():
        rel = elevation - sea_level
        high_cut = float(np.percentile(rel[valid], 58))
    else:
        rel = elevation - sea_level
        high_cut = 0.0
    highland_aligned = divide_cells & (
        (rel >= high_cut)
        | np.isin(province, [1, 5, 6, 7, 9])
    )
    divide_count = int(np.count_nonzero(divide_cells))
    divide_alignment_fraction = (
        float(np.count_nonzero(boundary_aligned | highland_aligned) / divide_count)
        if divide_count else 0.0
    )
    highland_alignment_fraction = (
        float(np.count_nonzero(highland_aligned) / divide_count)
        if divide_count else 0.0
    )

    path_crossings = 0
    downhill_violations = 0
    path_checks = 0
    for start in np.where(valid)[0]:
        c = int(start)
        bid = int(basin[c])
        previous_elev = float(flow_surface[c])
        for _ in range(96):
            nxt = int(receiver[c])
            if nxt == c:
                break
            if basin[nxt] >= 0 and int(basin[nxt]) != bid:
                path_crossings += 1
                break
            if flow_surface[nxt] > previous_elev + 1.0:
                downhill_violations += 1
                break
            previous_elev = float(flow_surface[nxt])
            c = nxt
        path_checks += 1
    flow_consistency_fraction = (
        float((path_checks - path_crossings) / path_checks) if path_checks else 0.0
    )
    downhill_fraction = (
        float((path_checks - downhill_violations) / path_checks) if path_checks else 0.0
    )
    drainage_basin_objects = list(objects.get("terrain.drainage_basins", []))
    drainage_divide_objects = list(objects.get("terrain.drainage_divides", []))
    acceptance = {
        "current_elevation_field_available": "terrain.elevation_m" in fields,
        "current_landform_context_present": (
            all(count > 0 for count in context_groups.values())
            or all(count > 0 for count in production_context_groups.values())
        ),
        "current_expected_residuals_recorded": expected_residuals_recorded,
        "production_drainage_fields_available": not missing_items,
        "production_drainage_objects_available": (
            len(drainage_basin_objects) > 0 and len(drainage_divide_objects) > 0
        ),
        "production_basins_contiguous": (
            major_component_count > 0 and component_failures == 0
        ),
        "production_divides_aligned": (
            divide_count > 0 and divide_alignment_fraction >= 0.70
        ),
        "production_flow_paths_stay_in_basins": flow_consistency_fraction >= 0.98,
        "production_flow_paths_mostly_downhill": downhill_fraction >= 0.98,
    }
    metrics = {
        "missing_item_count": int(len(missing_items)),
        "landform_object_count": int(len(landforms)),
        "landform_kind_count": int(len(kind_counts)),
        "divide_parent_context_count": int(context_groups["divide_parent_context"]),
        "basin_sink_context_count": int(context_groups["basin_sink_context"]),
        "margin_sink_context_count": int(context_groups["margin_sink_context"]),
        "production_divide_parent_context_count": int(
            production_context_groups["divide_parent_context"]),
        "production_basin_sink_context_count": int(
            production_context_groups["basin_sink_context"]),
        "production_margin_sink_context_count": int(
            production_context_groups["margin_sink_context"]),
        "drainage_basin_field_count": int(len(basin_ids)),
        "drainage_basin_object_count": int(len(drainage_basin_objects)),
        "drainage_divide_object_count": int(len(drainage_divide_objects)),
        "divide_cell_count": int(divide_count),
        "divide_fraction_of_land": (
            float(area[divide_cells].sum() / max(float(area[valid].sum()), 1.0))
            if valid.any() else 0.0
        ),
        "current_divide_alignment_fraction": float(divide_alignment_fraction),
        "current_highland_alignment_fraction": float(highland_alignment_fraction),
        "major_basin_component_count": int(major_component_count),
        "major_basin_component_failure_count": int(component_failures),
        "flow_path_check_count": int(path_checks),
        "flow_path_crossing_count": int(path_crossings),
        "flow_path_downhill_violation_count": int(downhill_violations),
        "current_flow_to_sink_consistency_fraction": float(flow_consistency_fraction),
        "current_downhill_path_fraction": float(downhill_fraction),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_drainage_context_audit_ready"
            if all(acceptance.values())
            else "generated_world_drainage_context_audit_incomplete"
        ),
        "missing_drainage_items": tuple(missing_items),
        "expected_current_residual_items": EXPECTED_CURRENT_RESIDUAL_ITEMS,
        "kind_counts": dict(sorted((kind, int(count)) for kind, count in kind_counts.items())),
        "limitations": {
            "production_drainage_objects_missing": bool(missing_items),
        },
        "metrics": metrics,
        "acceptance": acceptance,
    }


def drainage_divide_province_alignment_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_drainage_divide_summary()
    current = (
        current_generated_drainage_audit(world)
        if world is not None
        else {
            "status": "generated_world_drainage_context_audit_not_run",
            "acceptance": {
                "current_elevation_field_available": False,
                "current_landform_context_present": False,
                "current_expected_residuals_recorded": False,
                "production_drainage_fields_available": False,
                "production_drainage_objects_available": False,
                "production_basins_contiguous": False,
                "production_divides_aligned": False,
                "production_flow_paths_stay_in_basins": False,
                "production_flow_paths_mostly_downhill": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_divide_fixture_ready": reference["status"]
        == "drainage_divide_reference_ready",
        "divide_boundary_alignment": reference["acceptance"][
            "divide_boundary_alignment"],
        "flow_to_expected_sinks": reference["acceptance"]["flow_to_expected_sinks"],
        "basins_not_checkerboarded": reference["acceptance"][
            "basins_not_checkerboarded"],
        "current_generated_audit_available": current["status"]
        == "generated_world_drainage_context_audit_ready",
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_drainage_fields_available": current["acceptance"][
            "production_drainage_fields_available"],
        "production_drainage_objects_available": current["acceptance"][
            "production_drainage_objects_available"],
        "production_basins_contiguous": current["acceptance"][
            "production_basins_contiguous"],
        "production_divides_aligned": current["acceptance"][
            "production_divides_aligned"],
        "production_flow_paths_stay_in_basins": current["acceptance"][
            "production_flow_paths_stay_in_basins"],
    }
    return {
        "schema": SCHEMA,
        "status": (
            "drainage_divide_province_alignment_ready"
            if all(acceptance.values())
            else "drainage_divide_province_alignment_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
