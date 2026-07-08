"""P173.1 frame-level unsupported open-ocean shoal attribution."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.diagnostics.historical_geomorphology import (
    _dilate_mask,
    _distance_from_sources,
    _field,
    _frame_object_mask,
    _neighbor_range,
)
from aevum.modules.tectonics import (
    DOMAIN_ACCRETED_TERRANE,
    DOMAIN_LIP,
    ORIGIN_PLUME_IMPACT,
)
from aevum.modules.terrain import (
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RESTRICTED,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_RISE,
    OCEAN_DEPTH_SHELF,
    OCEAN_DEPTH_SLOPE,
    OCEAN_DEPTH_TRENCH,
    OCEAN_MARGIN_ACTIVE,
    OCEAN_MARGIN_OPEN,
    OCEAN_MARGIN_PASSIVE,
    OCEAN_MARGIN_RESTRICTED,
    OCEAN_MARGIN_RIDGE,
)


SCHEMA = "aevum.p173_unsupported_shoal_attribution.v1"

FRAME_METRIC_KEYS = (
    "cleanup_candidate_fraction_of_ocean",
    "structural_preserve_fraction_of_candidate",
    "object_support_fraction_of_candidate",
    "semantic_support_fraction_of_candidate",
    "post_cleanup_residual_fraction_of_ocean",
    "post_cleanup_residual_fraction_of_candidate",
)

P1732_FRAME_GLOBAL_KEYS = (
    "terrain.last_p1732_young_open_ocean_age_depth_used",
    "terrain.last_p1732_young_open_ocean_age_depth_land_mask_preserved",
    "terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction",
    "terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction",
    "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_before_m",
    "terrain.last_p1732_young_open_ocean_age_depth_mean_depth_after_m",
)

DEPTH_PROVINCE_NAMES = {
    0: "land",
    int(OCEAN_DEPTH_SHELF): "shelf",
    int(OCEAN_DEPTH_SLOPE): "slope",
    int(OCEAN_DEPTH_RISE): "rise",
    int(OCEAN_DEPTH_ABYSS): "abyss",
    int(OCEAN_DEPTH_RIDGE): "ridge",
    int(OCEAN_DEPTH_TRENCH): "trench",
    int(OCEAN_DEPTH_RESTRICTED): "restricted",
}

MARGIN_TYPE_NAMES = {
    0: "none",
    int(OCEAN_MARGIN_PASSIVE): "passive",
    int(OCEAN_MARGIN_ACTIVE): "active",
    int(OCEAN_MARGIN_RIDGE): "ridge",
    int(OCEAN_MARGIN_RESTRICTED): "restricted",
    int(OCEAN_MARGIN_OPEN): "open",
}


def unsupported_shoal_attribution_summary(
    world: Any,
    archive: Any,
    *,
    shallow_depth_m: float = 1500.0,
    coast_distance_steps: int = 8,
) -> dict[str, Any]:
    """Summarize per-frame cleanup ownership for unsupported shoals."""
    grid = getattr(world, "grid", None) or getattr(
        getattr(archive, "world", None), "grid", None)
    frames = list(getattr(archive, "frames", []) or [])
    if grid is None:
        return {
            "schema": SCHEMA,
            "frame_count": int(len(frames)),
            "usable_frame_count": 0,
            "skip_reason": "archive_grid_not_available",
            "frame_rows": [],
            "metric_extremes": _metric_extremes([]),
            "p1732_metric_extremes": _p1732_metric_extremes([]),
            "peak_residual_frame": {},
            "acceptance": {
                "attribution_completed": False,
                "generation_behavior_changed": False,
            },
        }

    rows = [
        frame_unsupported_shoal_attribution(
            grid,
            frame,
            shallow_depth_m=float(shallow_depth_m),
            coast_distance_steps=int(coast_distance_steps),
        )
        for frame in frames
    ]
    usable = [row for row in rows if bool(row.get("usable", False))]
    peak = (
        max(
            usable,
            key=lambda row: float(
                row.get("post_cleanup_residual_fraction_of_ocean", 0.0)),
        )
        if usable else {}
    )
    return {
        "schema": SCHEMA,
        "frame_count": int(len(rows)),
        "usable_frame_count": int(len(usable)),
        "config": {
            "shallow_depth_m": float(shallow_depth_m),
            "coast_distance_steps": int(coast_distance_steps),
        },
        "frame_rows": rows,
        "metric_extremes": _metric_extremes(usable),
        "p1732_metric_extremes": _p1732_metric_extremes(usable),
        "peak_residual_frame": _compact_peak_frame(peak),
        "acceptance": {
            "attribution_completed": bool(usable),
            "generation_behavior_changed": False,
        },
    }


def write_p173_unsupported_shoal_attribution(
    world: Any,
    archive: Any,
    outdir: str | Path,
    *,
    filename: str = "p173_unsupported_shoal_attribution.json",
) -> dict[str, Any]:
    """Write the P173.1 attribution JSON and return the same summary."""
    summary = unsupported_shoal_attribution_summary(world, archive)
    path = Path(outdir)
    path.mkdir(parents=True, exist_ok=True)
    (path / filename).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def frame_unsupported_shoal_attribution(
    grid: Any,
    frame: Any,
    *,
    shallow_depth_m: float = 1500.0,
    coast_distance_steps: int = 8,
) -> dict[str, Any]:
    """Return one frame's candidate/preserve/support/residual ownership."""
    fields = dict(getattr(frame, "fields", {}) or {})
    globals_ = dict(getattr(frame, "globals", {}) or {})
    n = int(getattr(grid, "n", 0))
    area = np.asarray(getattr(grid, "cell_area", np.ones(n)), dtype=np.float64)
    missing: list[str] = []
    elev, elev_ok = _field(fields, "terrain.elevation_m", n, missing)
    if not elev_ok:
        return {
            "schema": "aevum.p173_unsupported_shoal_frame.v1",
            "time_myr": float(getattr(frame, "time_myr", 0.0)),
            "usable": False,
            "skip_reason": "terrain.elevation_m_missing_or_invalid",
            "missing_fields": missing,
        }

    sea_level = float(globals_.get("ocean.sea_level_m", 0.0))
    rel = elev - sea_level
    land = rel >= 0.0
    ocean = ~land
    ocean_area = max(float(area[ocean].sum()), 1.0e-12)
    total_area = max(float(area.sum()), 1.0e-12)
    depth = np.maximum(-rel, 0.0)
    depth_province, _ = _field(fields, "ocean.depth_province", n, missing, default=0.0)
    margin_type, _ = _field(fields, "ocean.margin_type", n, missing, default=0.0)
    shelf_width, _ = _field(fields, "ocean.shelf_width", n, missing, default=-1.0)
    crust_age, has_age = _field(fields, "crust.age_myr", n, missing, default=np.nan)
    crust_type, _ = _field(fields, "crust.type", n, missing, default=0.0)
    crust_domain, _ = _field(fields, "crust.domain", n, missing, default=0.0)
    crust_origin, _ = _field(fields, "crust.origin", n, missing, default=0.0)
    rift_stage, _ = _field(fields, "terrain.rift_margin_stage", n, missing, default=0.0)
    if not ocean.any():
        return _empty_frame_row(frame, total_area=total_area)

    coast_distance = _distance_from_sources(
        grid, land_mask=land, domain=ocean, max_passes=int(coast_distance_steps))
    far_ocean = ocean & ((coast_distance < 0) | (coast_distance >= 4))
    candidate = far_ocean & (depth < float(shallow_depth_m))

    masks = _ownership_masks(
        grid,
        frame,
        ocean=ocean,
        depth=depth,
        depth_province=depth_province.astype(int),
        margin_type=margin_type.astype(int),
        coast_distance=coast_distance,
        crust_age=crust_age,
        has_age=has_age,
        crust_type=crust_type,
        crust_domain=crust_domain,
        crust_origin=crust_origin,
    )
    structural_preserve = masks["structural_preserve"]
    object_support = masks["object_support"]
    semantic_support = masks["semantic_support"]
    support = object_support | semantic_support
    post_cleanup_residual = candidate & ~structural_preserve & ~support

    candidate_area = max(float(area[candidate].sum()), 1.0e-12)
    residual_area = float(area[post_cleanup_residual].sum())
    return {
        "schema": "aevum.p173_unsupported_shoal_frame.v1",
        "time_myr": float(getattr(frame, "time_myr", 0.0)),
        "usable": True,
        "missing_fields": sorted(set(missing)),
        "ocean_area_fraction_of_world": float(area[ocean].sum() / total_area),
        "cleanup_candidate_fraction_of_ocean": float(area[candidate].sum() / ocean_area),
        "structural_preserve_fraction_of_candidate": _share(
            area, candidate & structural_preserve, candidate_area),
        "object_support_fraction_of_candidate": _share(
            area, candidate & object_support, candidate_area),
        "semantic_support_fraction_of_candidate": _share(
            area, candidate & semantic_support, candidate_area),
        "post_cleanup_residual_fraction_of_ocean": float(residual_area / ocean_area),
        "post_cleanup_residual_fraction_of_candidate": float(residual_area / candidate_area),
        "mask_counts": {
            "cleanup_candidate": int(np.count_nonzero(candidate)),
            "structural_preserve": int(np.count_nonzero(candidate & structural_preserve)),
            "object_support": int(np.count_nonzero(candidate & object_support)),
            "semantic_support": int(np.count_nonzero(candidate & semantic_support)),
            "post_cleanup_residual": int(np.count_nonzero(post_cleanup_residual)),
        },
        "mask_area_fraction_of_ocean": {
            "cleanup_candidate": float(area[candidate].sum() / ocean_area),
            "structural_preserve": float(area[candidate & structural_preserve].sum() / ocean_area),
            "object_support": float(area[candidate & object_support].sum() / ocean_area),
            "semantic_support": float(area[candidate & semantic_support].sum() / ocean_area),
            "post_cleanup_residual": float(residual_area / ocean_area),
        },
        "mask_fingerprints": {
            "cleanup_candidate": _mask_digest(candidate),
            "structural_preserve": _mask_digest(candidate & structural_preserve),
            "object_support": _mask_digest(candidate & object_support),
            "semantic_support": _mask_digest(candidate & semantic_support),
            "post_cleanup_residual": _mask_digest(post_cleanup_residual),
        },
        "p1732_young_open_ocean_depth_floor": _frame_global_metrics(
            globals_, P1732_FRAME_GLOBAL_KEYS),
        "residual_attribution": {
            "owner_hint": _owner_hint(
                post_cleanup_residual,
                area,
                depth,
                crust_age,
                depth_province.astype(int),
                margin_type.astype(int),
                rift_stage.astype(int),
            ),
            "depth_percentiles_m": _percentiles(depth, post_cleanup_residual),
            "crust_age_percentiles_myr": _percentiles(crust_age, post_cleanup_residual),
            "by_depth_province": _category_area_rows(
                depth_province.astype(int),
                post_cleanup_residual,
                area,
                DEPTH_PROVINCE_NAMES,
            ),
            "by_margin_type": _category_area_rows(
                margin_type.astype(int),
                post_cleanup_residual,
                area,
                MARGIN_TYPE_NAMES,
            ),
            "by_shelf_width": _category_area_rows(
                shelf_width.astype(int),
                post_cleanup_residual,
                area,
                {},
                max_rows=10,
            ),
            "by_rift_margin_stage": _category_area_rows(
                rift_stage.astype(int),
                post_cleanup_residual,
                area,
                {},
                max_rows=10,
            ),
            "by_coast_distance": _category_area_rows(
                coast_distance.astype(int),
                post_cleanup_residual,
                area,
                {},
                max_rows=10,
            ),
            "component_summary": _component_summary(
                grid,
                post_cleanup_residual,
                area,
                depth,
                crust_age,
            ),
        },
    }


def _ownership_masks(
    grid: Any,
    frame: Any,
    *,
    ocean: np.ndarray,
    depth: np.ndarray,
    depth_province: np.ndarray,
    margin_type: np.ndarray,
    coast_distance: np.ndarray,
    crust_age: np.ndarray,
    has_age: bool,
    crust_type: np.ndarray,
    crust_domain: np.ndarray,
    crust_origin: np.ndarray,
) -> dict[str, np.ndarray]:
    n = int(getattr(grid, "n", 0))
    object_ridge = _frame_object_mask(
        frame, "terrain.ocean_fabric",
        {"spreading_center", "ridge_segment"}, n) & ocean
    object_fracture = _frame_object_mask(
        frame, "terrain.ocean_fabric",
        {"transform_fault", "fracture_zone"}, n) & ocean
    object_trench = _frame_object_mask(
        frame, "terrain.margin_landforms",
        {"trench"}, n) & ocean
    object_plateau = _frame_object_mask(
        frame, "terrain.arc_plume_landforms",
        {"oceanic_plateau", "large_igneous_province"}, n) & ocean
    object_microcontinent = _frame_object_mask(
        frame, "terrain.arc_plume_landforms",
        {"microcontinent", "accreted_terrane"}, n) & ocean
    object_parented_shoal = _frame_object_mask(
        frame, "terrain.arc_plume_landforms",
        {"seamount_chain", "hotspot_track", "island_arc", "back_arc_basin"}, n) & ocean
    object_margin_parented_shoal = _frame_object_mask(
        frame,
        "terrain.margin_landforms",
        {"volcanic_arc", "forearc_accretionary_prism", "passive_margin_wedge", "delta_fan"},
        n,
    ) & ocean
    object_rift_margin = _frame_object_mask(
        frame, "terrain.rift_margin_sequences", {"rift_margin_sequence"}, n) & ocean
    object_support = _dilate_mask(
        grid,
        (
            object_ridge
            | object_fracture
            | object_trench
            | object_plateau
            | object_microcontinent
            | object_parented_shoal
            | object_margin_parented_shoal
            | object_rift_margin
        ),
        passes=1,
    ) & ocean

    ridge = ocean & (
        (depth_province == int(OCEAN_DEPTH_RIDGE))
        | (
            has_age
            & np.isfinite(crust_age)
            & (crust_age <= 15.0)
            & (depth <= 3300.0)
        )
        | object_ridge
    )
    trench = ocean & (
        (depth_province == int(OCEAN_DEPTH_TRENCH))
        | object_trench
    )
    age_contrast = _neighbor_range(grid, crust_age) if has_age else np.zeros(n, dtype=np.float64)
    if has_age and np.any(np.isfinite(age_contrast[ocean])):
        valid = age_contrast[ocean & np.isfinite(age_contrast)]
        cut = max(25.0, float(np.percentile(valid, 90))) if valid.size else 25.0
        fracture = (
            ocean
            & ~ridge
            & ~trench
            & ((coast_distance < 0) | (coast_distance >= 4))
            & np.isfinite(crust_age)
            & (crust_age >= 18.0)
            & (age_contrast >= cut)
        )
    else:
        fracture = np.zeros(n, dtype=bool)
    fracture |= object_fracture
    restricted = ocean & (
        (depth_province == int(OCEAN_DEPTH_RESTRICTED))
        | (margin_type == int(OCEAN_MARGIN_RESTRICTED))
    )
    structural_preserve = ridge | trench | fracture | restricted
    semantic_plateau = ocean & (
        (crust_domain.astype(int) == int(DOMAIN_LIP))
        | (crust_origin.astype(int) == int(ORIGIN_PLUME_IMPACT))
    )
    semantic_microcontinent = (
        ocean
        & (
            (crust_type >= 0.5)
            | (crust_domain.astype(int) == int(DOMAIN_ACCRETED_TERRANE))
        )
        & ((coast_distance < 0) | (coast_distance >= 3))
    )
    semantic_support = semantic_plateau | semantic_microcontinent
    return {
        "structural_preserve": structural_preserve,
        "object_support": object_support,
        "semantic_support": semantic_support,
    }


def _empty_frame_row(frame: Any, *, total_area: float) -> dict[str, Any]:
    return {
        "schema": "aevum.p173_unsupported_shoal_frame.v1",
        "time_myr": float(getattr(frame, "time_myr", 0.0)),
        "usable": True,
        "missing_fields": [],
        "ocean_area_fraction_of_world": 0.0,
        "cleanup_candidate_fraction_of_ocean": 0.0,
        "structural_preserve_fraction_of_candidate": 0.0,
        "object_support_fraction_of_candidate": 0.0,
        "semantic_support_fraction_of_candidate": 0.0,
        "post_cleanup_residual_fraction_of_ocean": 0.0,
        "post_cleanup_residual_fraction_of_candidate": 0.0,
        "mask_counts": {},
        "mask_area_fraction_of_ocean": {},
        "mask_fingerprints": {},
        "p1732_young_open_ocean_depth_floor": _frame_global_metrics(
            dict(getattr(frame, "globals", {}) or {}),
            P1732_FRAME_GLOBAL_KEYS,
        ),
        "residual_attribution": {
            "owner_hint": "no_ocean",
            "depth_percentiles_m": {},
            "crust_age_percentiles_myr": {},
            "by_depth_province": [],
            "by_margin_type": [],
            "by_shelf_width": [],
            "by_rift_margin_stage": [],
            "by_coast_distance": [],
            "component_summary": {"component_count": 0, "largest_components": []},
        },
    }


def _share(area: np.ndarray, mask: np.ndarray, denom_area: float) -> float:
    return float(area[np.asarray(mask, dtype=bool)].sum() / max(float(denom_area), 1.0e-12))


def _mask_digest(mask: np.ndarray) -> str:
    cells = np.where(np.asarray(mask, dtype=bool))[0].astype(np.int64)
    if cells.size == 0:
        return "empty"
    return hashlib.sha1(cells.tobytes()).hexdigest()[:16]


def _percentiles(values: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {}
    pct = np.percentile(vals, [0, 10, 25, 50, 75, 90, 100])
    labels = ("p0", "p10", "p25", "p50", "p75", "p90", "p100")
    return {label: float(value) for label, value in zip(labels, pct)}


def _category_area_rows(
    values: np.ndarray,
    mask: np.ndarray,
    area: np.ndarray,
    names: dict[int, str],
    *,
    max_rows: int = 12,
) -> list[dict[str, Any]]:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return []
    total = max(float(area[mask].sum()), 1.0e-12)
    rows: list[dict[str, Any]] = []
    ivals = np.asarray(values).astype(int)
    for value in np.unique(ivals[mask]):
        sub = mask & (ivals == int(value))
        rows.append({
            "value": int(value),
            "name": str(names.get(int(value), int(value))),
            "fraction_of_residual": float(area[sub].sum() / total),
            "cell_count": int(np.count_nonzero(sub)),
        })
    rows.sort(key=lambda row: (-float(row["fraction_of_residual"]), int(row["value"])))
    return rows[:int(max_rows)]


def _component_summary(
    grid: Any,
    mask: np.ndarray,
    area: np.ndarray,
    depth: np.ndarray,
    age: np.ndarray,
) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return {"component_count": 0, "largest_components": []}
    seen = np.zeros(mask.size, dtype=bool)
    total = max(float(area[mask].sum()), 1.0e-12)
    components: list[dict[str, Any]] = []
    for start in np.where(mask)[0]:
        if seen[int(start)]:
            continue
        stack = [int(start)]
        seen[int(start)] = True
        cells: list[int] = []
        while stack:
            cell = stack.pop()
            cells.append(cell)
            for nb in np.asarray(grid.neighbors[int(cell)], dtype=int):
                nb = int(nb)
                if mask[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        comp = np.asarray(cells, dtype=np.int64)
        comp_area = float(area[comp].sum())
        weights = area[comp]
        age_vals = age[comp]
        components.append({
            "cell_count": int(comp.size),
            "fraction_of_residual": float(comp_area / total),
            "centroid_lat": float(np.average(grid.lat[comp], weights=weights)),
            "centroid_lon": float(np.average(grid.lon[comp], weights=weights)),
            "depth_p50_m": float(np.median(depth[comp])),
            "crust_age_p50_myr": (
                float(np.nanmedian(age_vals)) if np.isfinite(age_vals).any() else 0.0
            ),
        })
    components.sort(key=lambda row: -float(row["fraction_of_residual"]))
    return {
        "component_count": int(len(components)),
        "largest_components": components[:8],
    }


def _owner_hint(
    residual: np.ndarray,
    area: np.ndarray,
    depth: np.ndarray,
    age: np.ndarray,
    depth_province: np.ndarray,
    margin_type: np.ndarray,
    rift_stage: np.ndarray,
) -> str:
    if not np.asarray(residual, dtype=bool).any():
        return "none"
    total = max(float(area[residual].sum()), 1.0e-12)

    def frac(mask: np.ndarray) -> float:
        return float(area[residual & mask].sum() / total)

    open_share = frac(margin_type == int(OCEAN_MARGIN_OPEN))
    abyss_rise_share = frac(
        (depth_province == int(OCEAN_DEPTH_ABYSS))
        | (depth_province == int(OCEAN_DEPTH_RISE))
    )
    rift_share = frac(rift_stage > 0)
    median_depth = float(np.median(depth[residual]))
    age_vals = age[residual]
    median_age = float(np.nanmedian(age_vals)) if np.isfinite(age_vals).any() else 0.0
    if open_share >= 0.55 and abyss_rise_share >= 0.55 and median_depth < 700.0:
        if 15.0 <= median_age <= 90.0:
            return "open_ocean_young_shallow_abyss_rise"
        return "open_ocean_shallow_abyss_rise"
    if rift_share >= 0.35:
        return "rift_margin_stage_shoal"
    if open_share >= 0.55:
        return "open_ocean_unattributed_shoal"
    return "mixed_unattributed_shoal"


def _metric_extremes(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for key in FRAME_METRIC_KEYS:
        values = np.asarray([
            float(row.get(key, 0.0)) for row in rows
        ], dtype=np.float64)
        out[key] = {
            "min": float(np.min(values)) if values.size else 0.0,
            "max": float(np.max(values)) if values.size else 0.0,
            "median": float(np.median(values)) if values.size else 0.0,
        }
    return out


def _p1732_metric_extremes(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for key in P1732_FRAME_GLOBAL_KEYS:
        values = np.asarray([
            float((row.get("p1732_young_open_ocean_depth_floor", {}) or {}).get(
                key, 0.0))
            for row in rows
        ], dtype=np.float64)
        out[key] = {
            "min": float(np.min(values)) if values.size else 0.0,
            "max": float(np.max(values)) if values.size else 0.0,
            "median": float(np.median(values)) if values.size else 0.0,
        }
    return out


def _frame_global_metrics(
    globals_: dict[str, Any],
    keys: tuple[str, ...],
) -> dict[str, float]:
    source = globals_ if isinstance(globals_, dict) else {}
    out: dict[str, float] = {}
    for key in keys:
        try:
            out[key] = float(source.get(key, 0.0))
        except (TypeError, ValueError):
            out[key] = 0.0
    return out


def _compact_peak_frame(row: dict[str, Any]) -> dict[str, Any]:
    if not row:
        return {}
    attribution = dict(row.get("residual_attribution", {}) or {})
    p1732 = dict(row.get("p1732_young_open_ocean_depth_floor", {}) or {})
    return {
        "time_myr": float(row.get("time_myr", 0.0)),
        "post_cleanup_residual_fraction_of_ocean": float(
            row.get("post_cleanup_residual_fraction_of_ocean", 0.0)),
        "post_cleanup_residual_fraction_of_candidate": float(
            row.get("post_cleanup_residual_fraction_of_candidate", 0.0)),
        "cleanup_candidate_fraction_of_ocean": float(
            row.get("cleanup_candidate_fraction_of_ocean", 0.0)),
        "owner_hint": str(attribution.get("owner_hint", "none")),
        "depth_percentiles_m": dict(attribution.get("depth_percentiles_m", {}) or {}),
        "crust_age_percentiles_myr": dict(
            attribution.get("crust_age_percentiles_myr", {}) or {}),
        "dominant_depth_province": _first_row(
            attribution.get("by_depth_province", [])),
        "dominant_margin_type": _first_row(attribution.get("by_margin_type", [])),
        "component_summary": dict(attribution.get("component_summary", {}) or {}),
        "p1732_young_open_ocean_depth_floor": p1732,
    }


def _first_row(rows: Any) -> dict[str, Any]:
    if isinstance(rows, list) and rows:
        return dict(rows[0])
    return {}
