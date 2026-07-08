"""Selected-snapshot high-resolution terrain refinement.

This module derives a high-resolution terminal snapshot from a P107 terminal
array archive.  It does not replay deep time.  Parent topology and process
semantics come from the source archive; the high-resolution layer adds a
process-conditioned elevation detail field for visual and diagnostic review.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aevum import render
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.p107_array_render import (
    _resolve_arrays_path,
    _resolve_metrics_path,
    render_p107_array_assets,
)
from aevum.modules.terrain import (
    CONT_DETAIL_ARC_MICROCONTINENT,
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    CONT_DETAIL_PLATEAU,
    CONT_DETAIL_RIFT_BASIN,
    CONT_DETAIL_SHIELD,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_RESTRICTED,
    OCEAN_DEPTH_TRENCH,
)


SCHEMA = "aevum.selected_snapshot_refinement.v1"


@dataclass(frozen=True)
class SelectedSnapshotRefinementConfig:
    source: str | Path
    outdir: str | Path
    target_cells: int = 72000
    width: int = 1600
    height: int = 800
    interpolation_k: int = 6
    detail_seed: int = 72000
    detail_strength: float = 1.0
    allow_process_islands: bool = False
    render_groups: tuple[str, ...] = ("all",)


_REFINEMENT_RENDER_GROUPS = frozenset({
    "all",
    "p107",
    "base",
    "hydrology",
    "marine",
    "shelf",
    "deep-ocean",
    "submarine",
    "island-atoll",
    "coastal",
})

_REFINEMENT_RENDER_GROUP_ALIASES = {
    "refinement": "base",
    "core": "base",
    "rivers": "hydrology",
    "river": "hydrology",
    "shelf-slope": "shelf",
    "shelf_slope": "shelf",
    "deep_ocean": "deep-ocean",
    "deep": "deep-ocean",
    "submarine-highlands": "submarine",
    "submarine_highlands": "submarine",
    "seamount": "submarine",
    "seamounts": "submarine",
    "island": "island-atoll",
    "islands": "island-atoll",
    "atoll": "island-atoll",
    "atolls": "island-atoll",
    "reef": "island-atoll",
    "reefs": "island-atoll",
    "reef-atoll": "island-atoll",
    "reef_atoll": "island-atoll",
    "coasts": "coastal",
    "coast": "coastal",
}


def _normalize_render_groups(groups: tuple[str, ...] | list[str] | None) -> frozenset[str]:
    if not groups:
        return frozenset({"all"})
    out: set[str] = set()
    for raw_group in groups:
        for raw_part in str(raw_group).split(","):
            part = raw_part.strip().lower().replace("_", "-")
            if not part:
                continue
            part = _REFINEMENT_RENDER_GROUP_ALIASES.get(part, part)
            if part not in _REFINEMENT_RENDER_GROUPS:
                valid = ", ".join(sorted(_REFINEMENT_RENDER_GROUPS))
                raise ValueError(f"unknown selected-snapshot render group {raw_part!r}; valid: {valid}")
            out.add(part)
    return frozenset(out or {"all"})


def _render_group_enabled(groups: frozenset[str], group: str) -> bool:
    if "all" in groups or group in groups:
        return True
    return "marine" in groups and group in {
        "shelf",
        "deep-ocean",
        "submarine",
        "island-atoll",
    }


def refine_selected_snapshot(
    config: SelectedSnapshotRefinementConfig,
) -> dict[str, Any]:
    """Build and render a high-resolution refinement package."""
    out = Path(config.outdir)
    out.mkdir(parents=True, exist_ok=True)
    render_groups = _normalize_render_groups(config.render_groups)

    metrics_path = _resolve_metrics_path(Path(config.source))
    parent_metrics = json.loads(metrics_path.read_text())
    arrays_path = _resolve_arrays_path(metrics_path, parent_metrics)
    with np.load(arrays_path) as loaded:
        source_arrays = {key: loaded[key] for key in loaded.files}

    source_n = int(source_arrays["grid_lat"].shape[0])
    target_n = int(config.target_cells)
    if target_n <= 0:
        raise ValueError("target_cells must be positive")
    source_grid = SphereGrid.fibonacci(source_n, CONSTANTS.EARTH_RADIUS)
    target_grid = SphereGrid.fibonacci(target_n, CONSTANTS.EARTH_RADIUS)
    nearest, interp_idx, interp_dist = _target_parent_indices(
        source_grid, target_grid, k=int(config.interpolation_k))

    sea_level = float(parent_metrics.get("sea_level_m", 0.0))
    fields = parent_metrics.get("array_archive", {}).get("manifest", {}).get(
        "fields", {})
    elevation_key = _field_key(fields, "terrain.elevation_m")
    parent_elev = np.asarray(source_arrays[elevation_key], dtype=np.float64)
    parent_rel_source = parent_elev - sea_level
    parent_rel = _interpolate_continuous(
        parent_rel_source,
        interp_idx,
        interp_dist,
        preserve_sign=True,
    )

    target_arrays = _resample_target_arrays(
        source_arrays,
        target_grid,
        nearest,
        interp_idx,
        interp_dist,
        fields,
        elevation_key,
    )

    context = _RefinementContext(
        target_grid=target_grid,
        arrays=target_arrays,
        fields=fields,
        parent_rel=parent_rel,
        sea_level=sea_level,
        detail_seed=int(config.detail_seed),
        detail_strength=float(config.detail_strength),
        allow_process_islands=bool(config.allow_process_islands),
    )
    refined_rel, detail_delta, amplitude = _refined_elevation(context)
    marine = _selected_snapshot_marine_microgeomorphology(context, refined_rel)
    refined_rel = marine.refined_rel
    detail_delta = detail_delta + marine.marine_delta
    micro = _selected_snapshot_microgeomorphology(context, refined_rel)
    refined_rel = micro.refined_rel
    detail_delta = detail_delta + micro.hydrology_delta
    coastal = _selected_snapshot_coastal_morphology(
        context,
        refined_rel,
        marine,
        micro,
    )
    refined_rel = coastal.refined_rel
    detail_delta = detail_delta + coastal.coastal_delta
    target_arrays[elevation_key] = (refined_rel + sea_level).astype(np.float32)
    target_arrays["field__selected_snapshot_parent_elevation_rel_m"] = (
        parent_rel.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_detail_delta_m"] = (
        detail_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_detail_amplitude_m"] = (
        amplitude.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_marine_delta_m"] = (
        marine.marine_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_ocean_coast_distance_passes"] = (
        marine.ocean_coast_distance.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_shelf_break_rank"] = (
        marine.shelf_break_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_shelf_slope_microrelief_delta_m"] = (
        marine.shelf_slope_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_deep_ocean_fabric_delta_m"] = (
        marine.deep_ocean_fabric_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_fracture_zone_rank"] = (
        marine.fracture_zone_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_abyssal_plain_fabric_rank"] = (
        marine.abyssal_plain_fabric_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_island_candidate_rank"] = (
        marine.island_candidate_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_atoll_candidate_rank"] = (
        marine.atoll_candidate_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_reef_rim_rank"] = (
        marine.reef_rim_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_atoll_lagoon_rank"] = (
        marine.atoll_lagoon_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_fringing_reef_rank"] = (
        marine.fringing_reef_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_process_island_promotion_rank"] = (
        marine.process_island_promotion_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_process_island_promotion_delta_m"] = (
        marine.process_island_promotion_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_islet_microshape_rank"] = (
        marine.islet_microshape_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_atoll_microshape_rank"] = (
        marine.atoll_microshape_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_island_atoll_microrelief_delta_m"] = (
        marine.island_atoll_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_submarine_highland_delta_m"] = (
        marine.submarine_highland_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_inland_seaway_tuning_delta_m"] = (
        np.asarray(marine.inland_seaway_tuning_delta, dtype=np.float64).astype(np.float32)
    )
    target_arrays["field__selected_snapshot_inland_seaway_tuning_rank"] = (
        np.asarray(marine.inland_seaway_tuning_rank, dtype=np.float64).astype(np.float32)
    )
    target_arrays["mask__selected_snapshot_inland_seaway_landback"] = (
        np.asarray(marine.inland_seaway_landback_mask, dtype=bool).astype(np.uint8)
    )
    target_arrays["field__selected_snapshot_seamount_peak_rank"] = (
        marine.seamount_peak_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_seamount_apron_rank"] = (
        marine.seamount_apron_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_oceanic_plateau_edge_rank"] = (
        marine.oceanic_plateau_edge_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_abyssal_hill_field_rank"] = (
        marine.abyssal_hill_field_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_hydrology_delta_m"] = (
        micro.hydrology_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_fluvial_microrelief_delta_m"] = (
        micro.fluvial_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_lowland_alluvial_microrelief_delta_m"] = (
        micro.lowland_alluvial_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_alluvial_fan_rank"] = (
        micro.alluvial_fan_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_lowland_plain_rank"] = (
        micro.lowland_plain_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_piedmont_apron_rank"] = (
        micro.piedmont_apron_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_delta_m"] = (
        coastal.coastal_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_process_microrelief_delta_m"] = (
        coastal.coastal_process_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_depositional_microrelief_delta_m"] = (
        coastal.coastal_depositional_microrelief_delta.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_plain_rank"] = (
        coastal.coastal_plain_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_cliff_rank"] = (
        coastal.coastal_cliff_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_shoreface_rank"] = (
        coastal.shoreface_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_barrier_lagoon_rank"] = (
        coastal.barrier_lagoon_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_estuary_rank"] = (
        coastal.estuary_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_delta_distributary_rank"] = (
        coastal.delta_distributary_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_estuary_funnel_rank"] = (
        coastal.estuary_funnel_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_barrier_spit_rank"] = (
        coastal.barrier_spit_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_delta_mouth_bar_rank"] = (
        coastal.delta_mouth_bar_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_estuary_tidal_channel_rank"] = (
        coastal.estuary_tidal_channel_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_coastal_depositional_plain_rank"] = (
        coastal.coastal_depositional_plain_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_strandplain_rank"] = (
        coastal.strandplain_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_tidal_flat_rank"] = (
        coastal.tidal_flat_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_flow_accumulation"] = (
        micro.flow_accumulation.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_river_rank"] = (
        micro.river_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_river_path_rank"] = (
        micro.river_path_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_basin_trunk_rank"] = (
        micro.basin_trunk_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_floodplain_rank"] = (
        micro.floodplain_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_meander_belt_rank"] = (
        micro.meander_belt_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_meander_scroll_rank"] = (
        micro.meander_scroll_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_floodplain_swale_rank"] = (
        micro.floodplain_swale_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_lake_basin_rank"] = (
        micro.lake_basin_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_lake_shoreline_rank"] = (
        micro.lake_shoreline_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_delta_fan_rank"] = (
        micro.delta_fan_rank.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_drainage_basin_id"] = (
        micro.drainage_basin_id.astype(np.int32)
    )
    target_arrays["field__selected_snapshot_land_coast_distance_passes"] = (
        micro.land_coast_distance.astype(np.float32)
    )
    target_arrays["field__selected_snapshot_river_receiver"] = (
        micro.river_receiver.astype(np.int32)
    )
    target_arrays["mask__selected_snapshot_lakes"] = micro.lake_mask.astype(np.uint8)
    target_arrays["mask__selected_snapshot_deltas"] = micro.delta_mask.astype(np.uint8)
    target_arrays["mask__selected_snapshot_delta_plain"] = (
        micro.delta_plain_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_reef_atoll"] = (
        marine.reef_atoll_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_marine_shoal"] = (
        marine.marine_shoal_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_seamount_shoal"] = (
        marine.seamount_shoal_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_oceanic_plateau_shoal"] = (
        marine.oceanic_plateau_shoal_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_microcontinent_shoal"] = (
        marine.microcontinent_shoal_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_island_arc_shoal"] = (
        marine.island_arc_shoal_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_inland_seaway_tuned"] = (
        np.asarray(marine.inland_seaway_tuning_mask, dtype=bool).astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_process_island_candidate"] = (
        marine.process_island_candidate_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_atoll_candidate"] = (
        marine.atoll_candidate_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_process_island_promoted"] = (
        marine.process_island_promoted_mask.astype(np.uint8)
    )
    target_arrays["mask__selected_snapshot_atoll_islet_promoted"] = (
        marine.atoll_islet_promoted_mask.astype(np.uint8)
    )

    arrays_out = out / "selected_snapshot_refined_arrays.npz"
    np.savez_compressed(arrays_out, **target_arrays)

    refined_metrics = _refined_metrics(
        config,
        metrics_path,
        arrays_path,
        arrays_out,
        parent_metrics,
        parent_rel,
        refined_rel,
        detail_delta,
        amplitude,
        marine,
        micro,
        coastal,
        context,
    )
    metrics_out = out / "selected_snapshot_refinement_metrics.json"
    metrics_out.write_text(
        json.dumps(refined_metrics, indent=2, default=_json_default) + "\n")

    p107_metrics = _p107_compatible_metrics(
        refined_metrics,
        parent_metrics,
        arrays_out,
        target_arrays,
        fields,
    )
    p107_metrics_out = out / "p107_terminal_metrics.json"
    p107_metrics_out.write_text(
        json.dumps(p107_metrics, indent=2, default=_json_default) + "\n")

    rendered_dir = out / "rendered"
    if _render_group_enabled(render_groups, "p107"):
        p107_summary = render_p107_array_assets(
            p107_metrics_out,
            rendered_dir,
            width=int(config.width),
            height=int(config.height),
        )
    else:
        rendered_dir.mkdir(parents=True, exist_ok=True)
        p107_summary = {"assets": {}}
    extra_assets = _render_refinement_assets(
        target_grid,
        parent_rel,
        refined_rel,
        detail_delta,
        amplitude,
        marine,
        micro,
        coastal,
        rendered_dir,
        width=int(config.width),
        height=int(config.height),
        render_groups=render_groups,
    )
    refined_metrics["assets"] = {
        "p107": p107_summary["assets"],
        "refinement": {key: str(path) for key, path in extra_assets.items()},
    }
    metrics_out.write_text(
        json.dumps(refined_metrics, indent=2, default=_json_default) + "\n")
    return refined_metrics


def render_selected_snapshot_refinement_assets(
    source: str | Path,
    *,
    render_groups: tuple[str, ...] = ("all",),
    width: int | None = None,
    height: int | None = None,
    outdir: str | Path | None = None,
) -> dict[str, Any]:
    """Render selected-snapshot QA assets from an existing refined array archive."""
    metrics_path = _resolve_selected_snapshot_metrics_path(Path(source))
    metrics = json.loads(metrics_path.read_text())
    source_dir = metrics_path.parent
    arrays_path = _resolve_selected_snapshot_arrays_path(metrics_path, metrics)
    with np.load(arrays_path) as loaded:
        arrays = {key: loaded[key] for key in loaded.files}

    n = int(arrays["grid_lat"].shape[0])
    grid = SphereGrid.fibonacci(n, CONSTANTS.EARTH_RADIUS)
    sea_level = float(metrics.get("sea_level_m", 0.0))
    width = int(width or metrics.get("width", 1600))
    height = int(height or metrics.get("height", 800))
    groups = _normalize_render_groups(render_groups)
    rendered_dir = Path(outdir) if outdir is not None else source_dir / "rendered"
    rendered_dir.mkdir(parents=True, exist_ok=True)

    refined_rel = _array_float(
        arrays, "field__terrain_elevation_m", n, default=0.0) - sea_level
    parent_rel = _array_float(
        arrays, "field__selected_snapshot_parent_elevation_rel_m", n, default=refined_rel)
    detail_delta = _array_float(
        arrays, "field__selected_snapshot_detail_delta_m", n, default=0.0)
    amplitude = _array_float(
        arrays, "field__selected_snapshot_detail_amplitude_m", n, default=0.0)
    marine = _marine_from_refined_arrays(arrays, refined_rel)
    micro = _micro_from_refined_arrays(arrays, refined_rel)
    coastal = _coastal_from_refined_arrays(arrays, refined_rel)

    p107_assets: dict[str, str] = {}
    p107_metrics_path = source_dir / "p107_terminal_metrics.json"
    if _render_group_enabled(groups, "p107") and p107_metrics_path.exists():
        p107_summary = render_p107_array_assets(
            p107_metrics_path,
            rendered_dir,
            width=width,
            height=height,
        )
        p107_assets = p107_summary.get("assets", {})

    refinement_assets = _render_refinement_assets(
        grid,
        parent_rel,
        refined_rel,
        detail_delta,
        amplitude,
        marine,
        micro,
        coastal,
        rendered_dir,
        width=width,
        height=height,
        render_groups=groups,
    )
    rendered_assets = {key: str(path) for key, path in refinement_assets.items()}
    if rendered_dir == source_dir / "rendered":
        existing_assets = metrics.get("assets", {})
        existing_refinement = dict(existing_assets.get("refinement", {}))
        existing_refinement.update(rendered_assets)
        existing_p107 = dict(existing_assets.get("p107", {}))
        existing_p107.update(p107_assets)
        metrics["assets"] = {
            "p107": existing_p107,
            "refinement": existing_refinement,
        }
        metrics_path.write_text(
            json.dumps(metrics, indent=2, default=_json_default) + "\n")

    summary = {
        "schema": f"{SCHEMA}.render_groups.v1",
        "source_metrics": str(metrics_path),
        "refined_arrays": str(arrays_path),
        "render_groups": sorted(groups),
        "width": width,
        "height": height,
        "assets": {
            "p107": p107_assets,
            "refinement": rendered_assets,
        },
    }
    summary_path = rendered_dir / "selected_snapshot_render_group_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, default=_json_default) + "\n")
    return summary


def _resolve_selected_snapshot_metrics_path(source: Path) -> Path:
    if source.is_dir():
        path = source / "selected_snapshot_refinement_metrics.json"
    else:
        path = source
    if not path.exists():
        raise FileNotFoundError(f"selected-snapshot metrics not found: {path}")
    return path


def _resolve_selected_snapshot_arrays_path(metrics_path: Path, metrics: dict[str, Any]) -> Path:
    raw = metrics.get("refined_arrays")
    if not raw:
        raise KeyError("selected-snapshot metrics missing refined_arrays")
    path = Path(str(raw))
    if path.exists():
        return path
    candidate = metrics_path.parent / path
    if candidate.exists():
        return candidate
    candidate = metrics_path.parent / path.name
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"selected-snapshot arrays not found: {raw}")


def _array_float(
    arrays: dict[str, np.ndarray],
    key: str,
    n: int,
    *,
    default: float | np.ndarray,
) -> np.ndarray:
    if key in arrays:
        return np.asarray(arrays[key], dtype=np.float64)
    if isinstance(default, np.ndarray):
        return np.asarray(default, dtype=np.float64).copy()
    return np.full(n, float(default), dtype=np.float64)


def _array_int(
    arrays: dict[str, np.ndarray],
    key: str,
    n: int,
    *,
    default: int,
) -> np.ndarray:
    if key in arrays:
        return np.asarray(arrays[key], dtype=np.int64)
    return np.full(n, int(default), dtype=np.int64)


def _array_bool(arrays: dict[str, np.ndarray], key: str, n: int) -> np.ndarray:
    if key in arrays:
        return np.asarray(arrays[key], dtype=bool)
    return np.zeros(n, dtype=bool)


def _marine_from_refined_arrays(
    arrays: dict[str, np.ndarray],
    refined_rel: np.ndarray,
) -> _MarineMicrogeomorphology:
    n = int(refined_rel.shape[0])
    return _MarineMicrogeomorphology(
        refined_rel=refined_rel,
        marine_delta=_array_float(arrays, "field__selected_snapshot_marine_delta_m", n, default=0.0),
        ocean_coast_distance=_array_float(
            arrays, "field__selected_snapshot_ocean_coast_distance_passes", n, default=-1.0),
        shelf_break_rank=_array_float(arrays, "field__selected_snapshot_shelf_break_rank", n, default=0.0),
        shelf_slope_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_shelf_slope_microrelief_delta_m", n, default=0.0),
        deep_ocean_fabric_delta=_array_float(
            arrays, "field__selected_snapshot_deep_ocean_fabric_delta_m", n, default=0.0),
        fracture_zone_rank=_array_float(arrays, "field__selected_snapshot_fracture_zone_rank", n, default=0.0),
        abyssal_plain_fabric_rank=_array_float(
            arrays, "field__selected_snapshot_abyssal_plain_fabric_rank", n, default=0.0),
        island_candidate_rank=_array_float(
            arrays, "field__selected_snapshot_island_candidate_rank", n, default=0.0),
        atoll_candidate_rank=_array_float(
            arrays, "field__selected_snapshot_atoll_candidate_rank", n, default=0.0),
        reef_atoll_mask=_array_bool(arrays, "mask__selected_snapshot_reef_atoll", n),
        marine_shoal_mask=_array_bool(arrays, "mask__selected_snapshot_marine_shoal", n),
        seamount_shoal_mask=_array_bool(arrays, "mask__selected_snapshot_seamount_shoal", n),
        oceanic_plateau_shoal_mask=_array_bool(
            arrays, "mask__selected_snapshot_oceanic_plateau_shoal", n),
        microcontinent_shoal_mask=_array_bool(
            arrays, "mask__selected_snapshot_microcontinent_shoal", n),
        island_arc_shoal_mask=_array_bool(arrays, "mask__selected_snapshot_island_arc_shoal", n),
        process_island_candidate_mask=_array_bool(
            arrays, "mask__selected_snapshot_process_island_candidate", n),
        atoll_candidate_mask=_array_bool(arrays, "mask__selected_snapshot_atoll_candidate", n),
        reef_rim_rank=_array_float(arrays, "field__selected_snapshot_reef_rim_rank", n, default=0.0),
        atoll_lagoon_rank=_array_float(arrays, "field__selected_snapshot_atoll_lagoon_rank", n, default=0.0),
        fringing_reef_rank=_array_float(arrays, "field__selected_snapshot_fringing_reef_rank", n, default=0.0),
        process_island_promotion_rank=_array_float(
            arrays, "field__selected_snapshot_process_island_promotion_rank", n, default=0.0),
        process_island_promotion_delta=_array_float(
            arrays, "field__selected_snapshot_process_island_promotion_delta_m", n, default=0.0),
        process_island_promoted_mask=_array_bool(
            arrays, "mask__selected_snapshot_process_island_promoted", n),
        atoll_islet_promoted_mask=_array_bool(
            arrays, "mask__selected_snapshot_atoll_islet_promoted", n),
        islet_microshape_rank=_array_float(
            arrays, "field__selected_snapshot_islet_microshape_rank", n, default=0.0),
        atoll_microshape_rank=_array_float(
            arrays, "field__selected_snapshot_atoll_microshape_rank", n, default=0.0),
        island_atoll_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_island_atoll_microrelief_delta_m", n, default=0.0),
        submarine_highland_delta=_array_float(
            arrays, "field__selected_snapshot_submarine_highland_delta_m", n, default=0.0),
        inland_seaway_tuning_delta=_array_float(
            arrays, "field__selected_snapshot_inland_seaway_tuning_delta_m", n, default=0.0),
        inland_seaway_tuning_mask=_array_bool(
            arrays, "mask__selected_snapshot_inland_seaway_tuned", n),
        inland_seaway_tuning_rank=_array_float(
            arrays, "field__selected_snapshot_inland_seaway_tuning_rank", n, default=0.0),
        inland_seaway_landback_mask=_array_bool(
            arrays, "mask__selected_snapshot_inland_seaway_landback", n),
        seamount_peak_rank=_array_float(arrays, "field__selected_snapshot_seamount_peak_rank", n, default=0.0),
        seamount_apron_rank=_array_float(arrays, "field__selected_snapshot_seamount_apron_rank", n, default=0.0),
        oceanic_plateau_edge_rank=_array_float(
            arrays, "field__selected_snapshot_oceanic_plateau_edge_rank", n, default=0.0),
        abyssal_hill_field_rank=_array_float(
            arrays, "field__selected_snapshot_abyssal_hill_field_rank", n, default=0.0),
    )


def _micro_from_refined_arrays(
    arrays: dict[str, np.ndarray],
    refined_rel: np.ndarray,
) -> _Microgeomorphology:
    n = int(refined_rel.shape[0])
    return _Microgeomorphology(
        refined_rel=refined_rel,
        hydrology_delta=_array_float(arrays, "field__selected_snapshot_hydrology_delta_m", n, default=0.0),
        fluvial_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_fluvial_microrelief_delta_m", n, default=0.0),
        lowland_alluvial_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_lowland_alluvial_microrelief_delta_m", n, default=0.0),
        flow_accumulation=_array_float(arrays, "field__selected_snapshot_flow_accumulation", n, default=0.0),
        river_rank=_array_float(arrays, "field__selected_snapshot_river_rank", n, default=0.0),
        river_path_rank=_array_float(arrays, "field__selected_snapshot_river_path_rank", n, default=0.0),
        basin_trunk_rank=_array_float(arrays, "field__selected_snapshot_basin_trunk_rank", n, default=0.0),
        floodplain_rank=_array_float(arrays, "field__selected_snapshot_floodplain_rank", n, default=0.0),
        meander_belt_rank=_array_float(arrays, "field__selected_snapshot_meander_belt_rank", n, default=0.0),
        meander_scroll_rank=_array_float(arrays, "field__selected_snapshot_meander_scroll_rank", n, default=0.0),
        floodplain_swale_rank=_array_float(arrays, "field__selected_snapshot_floodplain_swale_rank", n, default=0.0),
        alluvial_fan_rank=_array_float(arrays, "field__selected_snapshot_alluvial_fan_rank", n, default=0.0),
        lowland_plain_rank=_array_float(arrays, "field__selected_snapshot_lowland_plain_rank", n, default=0.0),
        piedmont_apron_rank=_array_float(arrays, "field__selected_snapshot_piedmont_apron_rank", n, default=0.0),
        lake_basin_rank=_array_float(arrays, "field__selected_snapshot_lake_basin_rank", n, default=0.0),
        lake_shoreline_rank=_array_float(arrays, "field__selected_snapshot_lake_shoreline_rank", n, default=0.0),
        delta_fan_rank=_array_float(arrays, "field__selected_snapshot_delta_fan_rank", n, default=0.0),
        drainage_basin_id=_array_int(arrays, "field__selected_snapshot_drainage_basin_id", n, default=-1),
        river_receiver=_array_int(arrays, "field__selected_snapshot_river_receiver", n, default=-1),
        lake_mask=_array_bool(arrays, "mask__selected_snapshot_lakes", n),
        delta_mask=_array_bool(arrays, "mask__selected_snapshot_deltas", n),
        delta_plain_mask=_array_bool(arrays, "mask__selected_snapshot_delta_plain", n),
        land_coast_distance=_array_float(
            arrays, "field__selected_snapshot_land_coast_distance_passes", n, default=-1.0),
    )


def _coastal_from_refined_arrays(
    arrays: dict[str, np.ndarray],
    refined_rel: np.ndarray,
) -> _CoastalMorphology:
    n = int(refined_rel.shape[0])
    return _CoastalMorphology(
        refined_rel=refined_rel,
        coastal_delta=_array_float(arrays, "field__selected_snapshot_coastal_delta_m", n, default=0.0),
        coastal_process_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_coastal_process_microrelief_delta_m", n, default=0.0),
        coastal_depositional_microrelief_delta=_array_float(
            arrays, "field__selected_snapshot_coastal_depositional_microrelief_delta_m", n, default=0.0),
        coastal_plain_rank=_array_float(arrays, "field__selected_snapshot_coastal_plain_rank", n, default=0.0),
        coastal_cliff_rank=_array_float(arrays, "field__selected_snapshot_coastal_cliff_rank", n, default=0.0),
        shoreface_rank=_array_float(arrays, "field__selected_snapshot_shoreface_rank", n, default=0.0),
        barrier_lagoon_rank=_array_float(arrays, "field__selected_snapshot_barrier_lagoon_rank", n, default=0.0),
        estuary_rank=_array_float(arrays, "field__selected_snapshot_estuary_rank", n, default=0.0),
        delta_distributary_rank=_array_float(
            arrays, "field__selected_snapshot_delta_distributary_rank", n, default=0.0),
        estuary_funnel_rank=_array_float(arrays, "field__selected_snapshot_estuary_funnel_rank", n, default=0.0),
        barrier_spit_rank=_array_float(arrays, "field__selected_snapshot_barrier_spit_rank", n, default=0.0),
        delta_mouth_bar_rank=_array_float(
            arrays, "field__selected_snapshot_delta_mouth_bar_rank", n, default=0.0),
        estuary_tidal_channel_rank=_array_float(
            arrays, "field__selected_snapshot_estuary_tidal_channel_rank", n, default=0.0),
        coastal_depositional_plain_rank=_array_float(
            arrays, "field__selected_snapshot_coastal_depositional_plain_rank", n, default=0.0),
        strandplain_rank=_array_float(arrays, "field__selected_snapshot_strandplain_rank", n, default=0.0),
        tidal_flat_rank=_array_float(arrays, "field__selected_snapshot_tidal_flat_rank", n, default=0.0),
    )


def _target_parent_indices(
    source_grid: SphereGrid,
    target_grid: SphereGrid,
    *,
    k: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kk = max(1, min(int(k), int(source_grid.n)))
    dist, idx = source_grid._kdtree.query(target_grid.xyz, k=kk)
    if kk == 1:
        idx = idx[:, None]
        dist = dist[:, None]
    return idx[:, 0].astype(np.int64), idx.astype(np.int64), dist.astype(np.float64)


def _resample_target_arrays(
    source_arrays: dict[str, np.ndarray],
    target_grid: SphereGrid,
    nearest: np.ndarray,
    interp_idx: np.ndarray,
    interp_dist: np.ndarray,
    fields: dict[str, str],
    elevation_key: str,
) -> dict[str, np.ndarray]:
    source_n = int(nearest.max()) + 1 if nearest.size else 0
    continuous_keys = {
        _field_key(fields, "crust.age_myr"),
        _field_key(fields, "ocean.shelf_width"),
        _field_key(fields, "tectonics.orogeny_age_myr"),
        _field_key(fields, "tectonics.volcanism_age_myr"),
    }
    target: dict[str, np.ndarray] = {
        "grid_lat": target_grid.lat.astype(np.float32),
        "grid_lon": target_grid.lon.astype(np.float32),
        "grid_cell_area_m2": target_grid.cell_area.astype(np.float64),
    }
    for key, arr in source_arrays.items():
        arr = np.asarray(arr)
        if key in {"grid_lat", "grid_lon", "grid_cell_area_m2"}:
            continue
        if arr.ndim != 1:
            continue
        if source_n and arr.shape[0] < source_n:
            continue
        if arr.shape[0] != np.asarray(source_arrays["grid_lat"]).shape[0]:
            continue
        if key == elevation_key:
            continue
        if key in continuous_keys and np.issubdtype(arr.dtype, np.floating):
            target[key] = _interpolate_continuous(
                arr.astype(np.float64),
                interp_idx,
                interp_dist,
            ).astype(arr.dtype)
        else:
            target[key] = arr[nearest].astype(arr.dtype, copy=False)
    return target


def _interpolate_continuous(
    values: np.ndarray,
    idx: np.ndarray,
    dist: np.ndarray,
    *,
    preserve_sign: bool = False,
) -> np.ndarray:
    vals = np.asarray(values, dtype=np.float64)[idx]
    nearest = vals[:, 0]
    finite = np.isfinite(vals)
    weights = 1.0 / np.maximum(dist, 1.0e-10) ** 2
    weights = np.where(finite, weights, 0.0)
    if preserve_sign:
        weights = np.where(np.signbit(vals) == np.signbit(nearest)[:, None],
                           weights, 0.0)
    weight_sum = weights.sum(axis=1)
    out = (vals * weights).sum(axis=1) / np.maximum(weight_sum, 1.0e-12)
    bad = weight_sum <= 0.0
    out[bad] = nearest[bad]
    if preserve_sign:
        land_flip = (nearest >= 0.0) & (out < 0.0)
        ocean_flip = (nearest < 0.0) & (out >= 0.0)
        out[land_flip | ocean_flip] = nearest[land_flip | ocean_flip]
    return out


@dataclass
class _RefinementContext:
    target_grid: SphereGrid
    arrays: dict[str, np.ndarray]
    fields: dict[str, str]
    parent_rel: np.ndarray
    sea_level: float
    detail_seed: int
    detail_strength: float
    allow_process_islands: bool


def _refined_elevation(
    context: _RefinementContext,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grid = context.target_grid
    rel = np.asarray(context.parent_rel, dtype=np.float64)
    land = rel >= 0.0
    ocean = ~land
    detail = _int_field(context, "terrain.continental_detail", default=0)
    detail_region = _int_field(
        context, "terrain.continental_detail_region_code", default=detail)
    ocean_depth = _int_field(context, "ocean.depth_province", default=0)
    crust_age = _float_field(context, "crust.age_myr", default=500.0)

    hierarchy = _int_field(context, "terrain.orogenic_parent_hierarchy", default=0)
    spine = _int_field(context, "terrain.orogenic_hierarchy_spine", default=0)
    halo = _int_field(context, "terrain.orogenic_shoulder_halo", default=0)
    apron = _int_field(context, "terrain.orogenic_highland_apron", default=0)

    ridge = _mask(context, "boundary__ridge") | _mask(
        context, "object__province_mid_ocean_ridge")
    trench = _mask(context, "boundary__trench") | _mask(
        context, "object__margin_trench")
    island_arc = _mask(context, "object__island_arc") | _mask(
        context, "object__volcanic_arc")
    seamount = _mask(context, "object__seamount_chain")
    plateau = _mask(context, "object__oceanic_plateau")
    microcontinent = _mask(context, "object__microcontinent")
    abyssal_hill = _mask(context, "object__abyssal_hill")

    smooth = _signed_texture(grid, context.detail_seed, scale=1)
    fine = _signed_texture(grid, context.detail_seed + 137, scale=2)
    ridge_texture = _ridged_texture(grid, context.detail_seed + 271)

    amplitude = np.zeros(grid.n, dtype=np.float64)
    delta = np.zeros(grid.n, dtype=np.float64)

    shield = detail_region == CONT_DETAIL_SHIELD
    platform = detail_region == CONT_DETAIL_PLATFORM
    basin = np.isin(detail_region, [CONT_DETAIL_BASIN, CONT_DETAIL_RIFT_BASIN])
    orogen = detail_region == CONT_DETAIL_OROGEN
    cont_plateau = detail_region == CONT_DETAIL_PLATEAU
    arc_micro = detail_region == CONT_DETAIL_ARC_MICROCONTINENT

    land_amp = np.zeros(grid.n, dtype=np.float64)
    land_amp += 28.0 * land
    land_amp += 35.0 * shield
    land_amp += 55.0 * platform
    land_amp += 45.0 * basin
    land_amp += 120.0 * orogen
    land_amp += 95.0 * cont_plateau
    land_amp += 85.0 * arc_micro
    land_amp += 95.0 * np.clip(hierarchy, 0, 3)
    land_amp += 145.0 * np.clip(spine, 0, 3)
    land_amp += 45.0 * (halo > 0)
    land_amp += 30.0 * (apron > 0)
    land_delta = land_amp * (0.62 * smooth + 0.38 * fine)
    land_delta += (hierarchy > 0) * (75.0 * np.clip(hierarchy, 0, 3) *
                                    (ridge_texture - 0.35))
    land_delta += (spine > 0) * (110.0 * np.clip(spine, 0, 3) *
                                 ridge_texture)
    land_delta -= basin * (35.0 + 35.0 * ridge_texture)
    lowland = land & (rel < 500.0) & ~orogen & (hierarchy <= 0)
    land_delta[lowland] *= 0.45

    young_ocean = np.exp(-np.clip(crust_age, 0.0, 250.0) / 75.0)
    ocean_amp = np.zeros(grid.n, dtype=np.float64)
    ocean_amp += 8.0 * ocean
    ocean_amp += 24.0 * young_ocean * ocean
    ocean_amp += 55.0 * abyssal_hill
    ocean_amp += 18.0 * (ocean_depth == OCEAN_DEPTH_ABYSS)
    ocean_amp += 135.0 * ridge
    ocean_amp += 220.0 * trench
    ocean_amp += 120.0 * seamount
    ocean_amp += 130.0 * plateau
    ocean_amp += 160.0 * island_arc
    ocean_amp += 130.0 * microcontinent
    ocean_delta = ocean_amp * (0.50 * smooth + 0.50 * fine)
    ocean_delta += ridge * (260.0 + 170.0 * ridge_texture)
    ocean_delta -= trench * (360.0 + 220.0 * ridge_texture)
    ocean_delta += (seamount | plateau | microcontinent | island_arc) * (
        130.0 + 280.0 * ridge_texture
    )
    ocean_delta += (ocean_depth == OCEAN_DEPTH_RIDGE) * 130.0
    ocean_delta -= (ocean_depth == OCEAN_DEPTH_TRENCH) * 260.0

    amplitude[land] = land_amp[land]
    amplitude[ocean] = ocean_amp[ocean]
    delta[land] = land_delta[land]
    delta[ocean] = ocean_delta[ocean]
    polar_damp = _polar_detail_damping(grid)
    amplitude *= polar_damp
    delta *= polar_damp
    delta *= float(context.detail_strength)

    refined = rel + delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean] = np.minimum(refined[ocean], -2.0)
    return refined, refined - rel, amplitude


@dataclass
class _MarineMicrogeomorphology:
    refined_rel: np.ndarray
    marine_delta: np.ndarray
    ocean_coast_distance: np.ndarray
    shelf_break_rank: np.ndarray
    shelf_slope_microrelief_delta: np.ndarray
    deep_ocean_fabric_delta: np.ndarray
    fracture_zone_rank: np.ndarray
    abyssal_plain_fabric_rank: np.ndarray
    island_candidate_rank: np.ndarray
    atoll_candidate_rank: np.ndarray
    reef_atoll_mask: np.ndarray
    marine_shoal_mask: np.ndarray
    seamount_shoal_mask: np.ndarray
    oceanic_plateau_shoal_mask: np.ndarray
    microcontinent_shoal_mask: np.ndarray
    island_arc_shoal_mask: np.ndarray
    process_island_candidate_mask: np.ndarray
    atoll_candidate_mask: np.ndarray
    reef_rim_rank: np.ndarray
    atoll_lagoon_rank: np.ndarray
    fringing_reef_rank: np.ndarray
    process_island_promotion_rank: np.ndarray
    process_island_promotion_delta: np.ndarray
    process_island_promoted_mask: np.ndarray
    atoll_islet_promoted_mask: np.ndarray
    islet_microshape_rank: np.ndarray
    atoll_microshape_rank: np.ndarray
    island_atoll_microrelief_delta: np.ndarray
    submarine_highland_delta: np.ndarray
    inland_seaway_tuning_delta: np.ndarray
    inland_seaway_tuning_mask: np.ndarray
    inland_seaway_tuning_rank: np.ndarray
    inland_seaway_landback_mask: np.ndarray
    seamount_peak_rank: np.ndarray
    seamount_apron_rank: np.ndarray
    oceanic_plateau_edge_rank: np.ndarray
    abyssal_hill_field_rank: np.ndarray


def _selected_snapshot_marine_microgeomorphology(
    context: _RefinementContext,
    refined_rel: np.ndarray,
) -> _MarineMicrogeomorphology:
    """Add conservative coastal shelf, reef, and seamount micro-bathymetry."""
    grid = context.target_grid
    rel = np.asarray(refined_rel, dtype=np.float64)
    land = rel >= 0.0
    ocean = ~land
    if not ocean.any():
        zeros = np.zeros(grid.n, dtype=np.float64)
        return _MarineMicrogeomorphology(
            refined_rel=rel.copy(),
            marine_delta=zeros,
            ocean_coast_distance=np.where(ocean, 0.0, -1.0),
            shelf_break_rank=zeros,
            shelf_slope_microrelief_delta=zeros,
            deep_ocean_fabric_delta=zeros,
            fracture_zone_rank=zeros,
            abyssal_plain_fabric_rank=zeros,
            island_candidate_rank=zeros,
            atoll_candidate_rank=zeros,
            reef_atoll_mask=np.zeros(grid.n, dtype=bool),
            marine_shoal_mask=np.zeros(grid.n, dtype=bool),
            seamount_shoal_mask=np.zeros(grid.n, dtype=bool),
            oceanic_plateau_shoal_mask=np.zeros(grid.n, dtype=bool),
            microcontinent_shoal_mask=np.zeros(grid.n, dtype=bool),
            island_arc_shoal_mask=np.zeros(grid.n, dtype=bool),
            process_island_candidate_mask=np.zeros(grid.n, dtype=bool),
            atoll_candidate_mask=np.zeros(grid.n, dtype=bool),
            reef_rim_rank=zeros,
            atoll_lagoon_rank=zeros,
            fringing_reef_rank=zeros,
            process_island_promotion_rank=zeros,
            process_island_promotion_delta=zeros,
            process_island_promoted_mask=np.zeros(grid.n, dtype=bool),
            atoll_islet_promoted_mask=np.zeros(grid.n, dtype=bool),
            islet_microshape_rank=zeros,
            atoll_microshape_rank=zeros,
            island_atoll_microrelief_delta=zeros,
            submarine_highland_delta=zeros,
            inland_seaway_tuning_delta=zeros,
            inland_seaway_tuning_mask=np.zeros(grid.n, dtype=bool),
            inland_seaway_tuning_rank=zeros,
            inland_seaway_landback_mask=np.zeros(grid.n, dtype=bool),
            seamount_peak_rank=zeros,
            seamount_apron_rank=zeros,
            oceanic_plateau_edge_rank=zeros,
            abyssal_hill_field_rank=zeros,
        )

    distance = _ocean_distance_from_land(grid, ocean)
    finite_ocean = ocean & np.isfinite(distance)
    ocean_depth = _int_field(context, "ocean.depth_province", default=0)
    crust_age = _float_field(context, "crust.age_myr", default=120.0)
    shelf_width = _float_field(context, "ocean.shelf_width", default=0.0)
    ridge = _mask(context, "boundary__ridge") | _mask(
        context, "object__province_mid_ocean_ridge")
    trench = _mask(context, "boundary__trench") | _mask(
        context, "object__margin_trench")
    active_margin = _mask(context, "boundary__active_margin") | _mask(
        context, "object__boundary_active_margin")
    passive_margin = _mask(context, "boundary__passive_margin") | _mask(
        context, "object__boundary_passive_margin")
    seamount = _mask(context, "object__seamount_chain")
    plateau = _mask(context, "object__oceanic_plateau")
    microcontinent = _mask(context, "object__microcontinent")
    abyssal_hill = _mask(context, "object__abyssal_hill")
    abyssal_plain = _mask(context, "object__abyssal_plain")
    transform = (
        _mask(context, "boundary__transform")
        | _mask(context, "object__boundary_transform")
        | _mask(context, "object__transform_fault")
    )
    fracture_zone = _mask(context, "object__fracture_zone") | transform
    island_arc = _mask(context, "object__island_arc") | _mask(
        context, "object__volcanic_arc")
    process_high = ocean & (seamount | plateau | microcontinent | island_arc)
    gateway_id = _int_field(context, "ocean.gateway_id", default=-1)
    gateway_system_id = _int_field(context, "ocean.gateway_system_id", default=-1)

    width = 1.45 + 0.34 * np.clip(shelf_width, 0.0, 10.0)
    width += np.where(passive_margin, 0.85, 0.0)
    width -= np.where(active_margin | trench, 0.55, 0.0)
    width = np.clip(width, 1.2, 5.4)
    d = np.where(finite_ocean, distance, 99.0)
    shelf_profile = -(80.0 + 95.0 * d + 56.0 * np.maximum(d - 1.0, 0.0) ** 1.35)
    shelf_domain = (
        finite_ocean
        & (d <= width)
        & (rel < -20.0)
        & (rel > -2400.0)
        & ~trench
        & (ocean_depth != OCEAN_DEPTH_TRENCH)
    )
    shelf_raise = np.zeros(grid.n, dtype=np.float64)
    too_deep_for_shelf = shelf_domain & (rel < shelf_profile)
    shelf_raise[too_deep_for_shelf] = np.minimum(
        360.0,
        0.22 * (shelf_profile[too_deep_for_shelf] - rel[too_deep_for_shelf]),
    )

    shelf_break_rank = np.zeros(grid.n, dtype=np.float64)
    shelf_break_rank[finite_ocean] = np.exp(-((d[finite_ocean] - width[finite_ocean])
                                              / 1.35) ** 2)
    shelf_break_rank *= np.where(trench | ridge, 0.35, 1.0)
    shelf_break_rank *= finite_ocean & (d <= width + 1.7) & (rel > -2800.0)
    shelf_break_rank *= ocean
    slope_domain = finite_ocean & (d > width) & (d <= width + 2.5) & ~ridge
    slope_deepen = np.zeros(grid.n, dtype=np.float64)
    slope_deepen[slope_domain] = -45.0 * shelf_break_rank[slope_domain] * np.clip(
        (rel[slope_domain] + 3600.0) / 2600.0,
        0.0,
        1.0,
    )

    ridge_texture = _ridged_texture(grid, context.detail_seed + 1409)
    fine = _signed_texture(grid, context.detail_seed + 1553, scale=4)
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    low_polar = lat_abs <= 72.0
    shoal_mask = process_high & low_polar & ~trench & (rel > -2400.0)
    shoal_delta = np.zeros(grid.n, dtype=np.float64)
    shoal_delta[shoal_mask] = 45.0 + 125.0 * ridge_texture[shoal_mask]
    shoal_delta[shoal_mask & plateau] += 35.0
    shoal_delta[shoal_mask & microcontinent] += 55.0
    shoal_delta[shoal_mask & island_arc] += 50.0

    tropical = lat_abs <= 34.0
    reef_domain = (
        ocean
        & tropical
        & process_high
        & _offshore_atoll_candidate_gate(
            d,
            seamount=seamount,
            plateau=plateau,
            microcontinent=microcontinent,
            island_arc=island_arc,
        )
        & (d <= 14.0)
        & (rel > -1800.0)
        & (rel < -8.0)
        & ~trench
    )
    reef_seed = np.zeros(grid.n, dtype=bool)
    if reef_domain.any():
        score = (
            1.00 * np.clip((rel + 1800.0) / 1800.0, 0.0, 1.0)
            + 0.48 * ridge_texture
            + 0.24 * np.clip((34.0 - lat_abs) / 34.0, 0.0, 1.0)
            + 0.18 * _open_ocean_candidate_weight(d, start=3.0, full=9.0)
            + 0.22 * microcontinent
            + 0.18 * island_arc
        )
        count = min(max(8, int(round(0.00055 * grid.n))),
                    int(np.count_nonzero(reef_domain)))
        chosen = np.flatnonzero(reef_domain)[
            np.argsort(score[reef_domain])[-count:]
        ]
        reef_seed[chosen] = True
    reef_atoll = reef_seed
    if reef_seed.any():
        reef_atoll = _dilate_mask(grid, reef_seed, passes=1) & reef_domain
        reef_atoll |= reef_seed
    reef_delta = np.zeros(grid.n, dtype=np.float64)
    reef_target = -28.0 - 12.0 * (1.0 - ridge_texture)
    reef_lift = reef_atoll & (rel < reef_target)
    reef_delta[reef_lift] = 0.62 * (reef_target[reef_lift] - rel[reef_lift])

    marine_delta = shelf_raise + slope_deepen + shoal_delta + reef_delta
    marine_delta += np.where(finite_ocean & (d <= 2.0), 16.0 * fine, 0.0)
    marine_delta[trench] = np.minimum(marine_delta[trench], 0.0)
    marine_delta[~ocean] = 0.0
    potential_rel = rel + marine_delta
    island_candidate_rank, process_island_candidate = _marine_island_candidates(
        grid,
        rel,
        potential_rel,
        process_high,
        trench,
        seamount,
        plateau,
        microcontinent,
        island_arc,
        ridge_texture,
        low_polar,
        d,
    )
    atoll_candidate_rank = np.zeros(grid.n, dtype=np.float64)
    if reef_atoll.any():
        atoll_base = (
            0.42
            + 0.30 * np.clip((potential_rel + 260.0) / 260.0, 0.0, 1.0)
            + 0.18 * ridge_texture
            + 0.10 * np.clip((34.0 - lat_abs) / 34.0, 0.0, 1.0)
        )
        atoll_candidate_rank[reef_atoll] = np.clip(atoll_base[reef_atoll], 0.0, 1.0)
    atoll_candidate = reef_atoll & (atoll_candidate_rank >= 0.48)
    reef_rim_rank, atoll_lagoon_rank, fringing_reef_rank = (
        _marine_reef_atoll_morphology(
            grid,
            rel,
            potential_rel,
            reef_domain,
            reef_atoll,
            atoll_candidate,
            atoll_candidate_rank,
            process_island_candidate,
            process_high,
            microcontinent,
            island_arc,
            plateau,
            trench,
            d,
            ridge_texture,
            tropical,
        )
    )
    reef_morph_delta = _reef_atoll_morphology_delta(
        rel,
        reef_rim_rank,
        atoll_lagoon_rank,
        fringing_reef_rank,
    )
    marine_delta += reef_morph_delta
    (
        submarine_highland_delta,
        seamount_peak_rank,
        seamount_apron_rank,
        oceanic_plateau_edge_rank,
        abyssal_hill_field_rank,
    ) = _marine_submarine_highland_morphology(
        grid,
        rel,
        rel + marine_delta,
        seamount,
        plateau,
        microcontinent,
        abyssal_hill,
        ridge,
        trench,
        ocean_depth,
        d,
        ridge_texture,
        fine,
        low_polar,
    )
    marine_delta += submarine_highland_delta
    process_support = (
        process_high
        | ridge
        | trench
        | active_margin
        | fracture_zone
    )
    (
        inland_seaway_tuning_delta,
        inland_seaway_tuning_mask,
        inland_seaway_tuning_rank,
        inland_seaway_landback_mask,
    ) = _unsupported_inland_shallow_seaway_tuning(
        grid,
        rel,
        ocean,
        d,
        ocean_depth,
        shelf_width,
        gateway_id,
        gateway_system_id,
        process_support,
    )
    marine_delta += inland_seaway_tuning_delta

    process_island_promotion_delta = np.zeros(grid.n, dtype=np.float64)
    process_island_promotion_rank = np.zeros(grid.n, dtype=np.float64)
    process_island_promoted = np.zeros(grid.n, dtype=bool)
    atoll_islet_promoted = np.zeros(grid.n, dtype=bool)
    if context.allow_process_islands:
        (
            process_island_promotion_delta,
            process_island_promotion_rank,
            process_island_promoted,
            atoll_islet_promoted,
        ) = _marine_process_island_promotion(
            grid,
            rel,
            rel + marine_delta,
            island_candidate_rank,
            atoll_candidate_rank,
            reef_rim_rank,
            atoll_lagoon_rank,
            fringing_reef_rank,
            process_island_candidate,
            atoll_candidate,
            seamount,
            plateau,
            microcontinent,
            island_arc,
            trench,
            low_polar,
        )
        marine_delta += process_island_promotion_delta

    promoted_land = process_island_promoted | inland_seaway_landback_mask
    refined = rel + marine_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean & ~promoted_land] = np.minimum(
        refined[ocean & ~promoted_land],
        -2.0,
    )
    refined[promoted_land] = np.maximum(
        refined[promoted_land],
        2.0,
    )
    if process_island_promoted.any():
        process_island_promotion_delta[process_island_promoted] = (
            refined[process_island_promoted] + 2.0
        )
    shelf_slope_baseline = refined.copy()
    shelf_slope_microrelief_delta = _marine_shelf_slope_microrelief_delta(
        grid,
        rel,
        refined,
        shelf_break_rank,
        d,
        passive_margin,
        active_margin,
        trench,
        ridge,
        promoted_land,
        low_polar,
        context.detail_seed,
    )
    marine_delta += shelf_slope_microrelief_delta
    refined = rel + marine_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean & ~promoted_land] = np.minimum(
        refined[ocean & ~promoted_land],
        -2.0,
    )
    refined[promoted_land] = np.maximum(
        refined[promoted_land],
        2.0,
    )
    shelf_slope_microrelief_delta = refined - shelf_slope_baseline
    deep_ocean_fabric_baseline = refined.copy()
    (
        deep_ocean_fabric_delta,
        fracture_zone_rank,
        abyssal_plain_fabric_rank,
    ) = _marine_deep_ocean_fabric_delta(
        grid,
        rel,
        refined,
        ocean_depth,
        crust_age,
        d,
        fracture_zone,
        transform,
        abyssal_plain,
        ridge,
        trench,
        promoted_land,
        low_polar,
        context.detail_seed,
    )
    marine_delta += deep_ocean_fabric_delta
    refined = rel + marine_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean & ~promoted_land] = np.minimum(
        refined[ocean & ~promoted_land],
        -2.0,
    )
    refined[promoted_land] = np.maximum(
        refined[promoted_land],
        2.0,
    )
    deep_ocean_fabric_delta = refined - deep_ocean_fabric_baseline
    islet_microshape_rank, atoll_microshape_rank = _marine_island_microshape_ranks(
        island_candidate_rank,
        atoll_candidate_rank,
        reef_rim_rank,
        atoll_lagoon_rank,
        fringing_reef_rank,
        process_island_candidate,
        atoll_candidate,
        process_island_promotion_rank,
        process_island_promoted,
        atoll_islet_promoted,
        ocean,
        low_polar,
    )
    baseline_refined = refined.copy()
    island_atoll_microrelief_delta = _marine_island_atoll_microrelief_delta(
        grid,
        rel,
        refined,
        islet_microshape_rank,
        atoll_microshape_rank,
        reef_rim_rank,
        atoll_lagoon_rank,
        fringing_reef_rank,
        process_island_promoted,
        atoll_islet_promoted,
        ocean,
        ridge_texture,
        low_polar,
        context.detail_seed,
    )
    marine_delta += island_atoll_microrelief_delta
    refined = rel + marine_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean & ~promoted_land] = np.minimum(
        refined[ocean & ~promoted_land],
        -2.0,
    )
    refined[promoted_land] = np.maximum(
        refined[promoted_land],
        2.0,
    )
    island_atoll_microrelief_delta = refined - baseline_refined
    (
        final_landback_delta,
        final_landback_mask,
    ) = _final_unsupported_shallow_corridor_landback(
        grid,
        rel,
        refined,
        ocean,
        process_support,
    )
    if np.any(final_landback_mask):
        marine_delta += final_landback_delta
        inland_seaway_tuning_delta += final_landback_delta
        inland_seaway_tuning_mask |= final_landback_mask
        inland_seaway_landback_mask |= final_landback_mask
        inland_seaway_tuning_rank[final_landback_mask] = np.maximum(
            inland_seaway_tuning_rank[final_landback_mask],
            1.0,
        )
        promoted_land |= final_landback_mask
        refined = rel + marine_delta
        refined[land] = np.maximum(refined[land], 2.0)
        refined[ocean & ~promoted_land] = np.minimum(
            refined[ocean & ~promoted_land],
            -2.0,
        )
        refined[promoted_land] = np.maximum(
            refined[promoted_land],
            2.0,
        )
    distance_out = np.where(finite_ocean, distance, -1.0)
    seamount_shoal = shoal_mask & seamount
    plateau_shoal = shoal_mask & plateau
    microcontinent_shoal = shoal_mask & microcontinent
    island_arc_shoal = shoal_mask & island_arc
    return _MarineMicrogeomorphology(
        refined_rel=refined,
        marine_delta=refined - rel,
        ocean_coast_distance=distance_out,
        shelf_break_rank=shelf_break_rank,
        shelf_slope_microrelief_delta=shelf_slope_microrelief_delta,
        deep_ocean_fabric_delta=deep_ocean_fabric_delta,
        fracture_zone_rank=fracture_zone_rank,
        abyssal_plain_fabric_rank=abyssal_plain_fabric_rank,
        island_candidate_rank=island_candidate_rank,
        atoll_candidate_rank=atoll_candidate_rank,
        reef_atoll_mask=reef_atoll,
        marine_shoal_mask=shoal_mask,
        seamount_shoal_mask=seamount_shoal,
        oceanic_plateau_shoal_mask=plateau_shoal,
        microcontinent_shoal_mask=microcontinent_shoal,
        island_arc_shoal_mask=island_arc_shoal,
        process_island_candidate_mask=process_island_candidate,
        atoll_candidate_mask=atoll_candidate,
        reef_rim_rank=reef_rim_rank,
        atoll_lagoon_rank=atoll_lagoon_rank,
        fringing_reef_rank=fringing_reef_rank,
        process_island_promotion_rank=process_island_promotion_rank,
        process_island_promotion_delta=process_island_promotion_delta,
        process_island_promoted_mask=process_island_promoted,
        atoll_islet_promoted_mask=atoll_islet_promoted,
        islet_microshape_rank=islet_microshape_rank,
        atoll_microshape_rank=atoll_microshape_rank,
        island_atoll_microrelief_delta=island_atoll_microrelief_delta,
        submarine_highland_delta=submarine_highland_delta,
        inland_seaway_tuning_delta=inland_seaway_tuning_delta,
        inland_seaway_tuning_mask=inland_seaway_tuning_mask,
        inland_seaway_tuning_rank=inland_seaway_tuning_rank,
        inland_seaway_landback_mask=inland_seaway_landback_mask,
        seamount_peak_rank=seamount_peak_rank,
        seamount_apron_rank=seamount_apron_rank,
        oceanic_plateau_edge_rank=oceanic_plateau_edge_rank,
        abyssal_hill_field_rank=abyssal_hill_field_rank,
    )


@dataclass
class _Microgeomorphology:
    refined_rel: np.ndarray
    hydrology_delta: np.ndarray
    fluvial_microrelief_delta: np.ndarray
    lowland_alluvial_microrelief_delta: np.ndarray
    flow_accumulation: np.ndarray
    river_rank: np.ndarray
    river_path_rank: np.ndarray
    basin_trunk_rank: np.ndarray
    floodplain_rank: np.ndarray
    meander_belt_rank: np.ndarray
    meander_scroll_rank: np.ndarray
    floodplain_swale_rank: np.ndarray
    alluvial_fan_rank: np.ndarray
    lowland_plain_rank: np.ndarray
    piedmont_apron_rank: np.ndarray
    lake_basin_rank: np.ndarray
    lake_shoreline_rank: np.ndarray
    delta_fan_rank: np.ndarray
    drainage_basin_id: np.ndarray
    river_receiver: np.ndarray
    lake_mask: np.ndarray
    delta_mask: np.ndarray
    delta_plain_mask: np.ndarray
    land_coast_distance: np.ndarray


def _selected_snapshot_microgeomorphology(
    context: _RefinementContext,
    refined_rel: np.ndarray,
) -> _Microgeomorphology:
    """Add selected-snapshot river, lake, and delta microgeomorphology.

    This is a lightweight terminal-snapshot proxy, not a deep-time hydrology
    model.  It uses inherited relief and process fields to make local drainage
    visible without changing the selected parent coastline.
    """
    grid = context.target_grid
    rel = np.asarray(refined_rel, dtype=np.float64)
    land = rel >= 0.0
    ocean = ~land
    distance = _land_distance_from_ocean(grid, land)
    finite_land = land & np.isfinite(distance)
    if not finite_land.any():
        zeros = np.zeros(grid.n, dtype=np.float64)
        return _Microgeomorphology(
            refined_rel=rel.copy(),
            hydrology_delta=zeros,
            fluvial_microrelief_delta=zeros,
            lowland_alluvial_microrelief_delta=zeros,
            flow_accumulation=zeros,
            river_rank=zeros,
            river_path_rank=zeros,
            basin_trunk_rank=zeros,
            floodplain_rank=zeros,
            meander_belt_rank=zeros,
            meander_scroll_rank=zeros,
            floodplain_swale_rank=zeros,
            alluvial_fan_rank=zeros,
            lowland_plain_rank=zeros,
            piedmont_apron_rank=zeros,
            lake_basin_rank=zeros,
            lake_shoreline_rank=zeros,
            delta_fan_rank=zeros,
            drainage_basin_id=np.zeros(grid.n, dtype=np.int32),
            river_receiver=np.full(grid.n, -1, dtype=np.int64),
            lake_mask=np.zeros(grid.n, dtype=bool),
            delta_mask=np.zeros(grid.n, dtype=bool),
            delta_plain_mask=np.zeros(grid.n, dtype=bool),
            land_coast_distance=np.where(land, 0.0, -1.0),
        )

    detail = _int_field(context, "terrain.continental_detail", default=0)
    detail_region = _int_field(
        context, "terrain.continental_detail_region_code", default=detail)
    hierarchy = _int_field(context, "terrain.orogenic_parent_hierarchy", default=0)
    shoulder = _int_field(context, "terrain.orogenic_shoulder_halo", default=0)
    apron = _int_field(context, "terrain.orogenic_highland_apron", default=0)
    basin = np.isin(detail_region, [CONT_DETAIL_BASIN, CONT_DETAIL_RIFT_BASIN])
    orogen = (detail_region == CONT_DETAIL_OROGEN) | (hierarchy > 0)
    coastal_land = finite_land & (distance <= 0.0)

    lat_r = np.radians(grid.lat)
    lat_factor = np.clip(np.cos(lat_r), 0.12, 1.0) ** 0.55
    coast_factor = 0.45 + 0.55 * np.exp(-np.maximum(distance, 0.0) / 26.0)
    high_relief_factor = 0.78 + 0.22 * np.clip(rel / 2500.0, 0.0, 1.0)
    runoff = np.where(
        finite_land,
        lat_factor * coast_factor * high_relief_factor,
        0.0,
    )
    runoff *= np.where(orogen, 1.20, 1.0)
    runoff *= np.where(basin, 0.82, 1.0)
    runoff *= grid.cell_area / max(float(np.mean(grid.cell_area)), 1.0)

    potential = rel + np.where(finite_land, distance, 0.0) * 72.0
    potential -= np.where(basin, 95.0, 0.0)
    potential += _signed_texture(grid, context.detail_seed + 901, scale=1) * 18.0
    potential[~finite_land] = np.inf

    receiver = _downslope_receivers(
        grid, finite_land, coastal_land, potential, distance)
    accumulation = runoff.copy()
    land_cells = np.flatnonzero(finite_land)
    order = land_cells[np.argsort(distance[land_cells])[::-1]]
    for cell in order:
        rec = int(receiver[int(cell)])
        if rec >= 0 and rec != int(cell):
            accumulation[rec] += accumulation[int(cell)]

    land_acc = accumulation[finite_land]
    river_rank = np.zeros(grid.n, dtype=np.float64)
    if land_acc.size:
        lo = float(np.nanpercentile(land_acc, 87.0))
        hi = max(float(np.nanpercentile(land_acc, 99.2)), lo + 1.0e-9)
        river_rank[finite_land] = np.clip((accumulation[finite_land] - lo) / (hi - lo),
                                          0.0, 1.0)
    river_rank[coastal_land & (accumulation > np.nanpercentile(land_acc, 82.0))] = np.maximum(
        river_rank[coastal_land & (accumulation > np.nanpercentile(land_acc, 82.0))],
        0.35,
    )
    river_rank = _connect_major_river_paths(
        receiver,
        finite_land,
        accumulation,
        river_rank,
        source_threshold=0.55,
    )
    river_path_rank = _major_river_path_rank(
        grid,
        receiver,
        finite_land,
        coastal_land,
        accumulation,
        distance,
        river_rank,
    )
    basin_id, basin_trunk_rank, floodplain_rank = _drainage_basin_objects(
        grid,
        receiver,
        finite_land,
        coastal_land,
        accumulation,
        distance,
        river_rank,
        rel,
        basin,
    )
    river_path_rank = np.maximum(river_path_rank, basin_trunk_rank)
    meander_belt_rank = _meander_belt_rank(
        grid,
        finite_land,
        rel,
        distance,
        basin_trunk_rank,
        floodplain_rank,
        basin,
    )

    sink = finite_land & (receiver == np.arange(grid.n))
    sink_acc_threshold = float(np.nanpercentile(land_acc, 88.0)) if land_acc.size else 0.0
    lake_candidates = (
        sink
        & (distance >= 5.0)
        & (accumulation >= sink_acc_threshold)
        & ((rel < 650.0) | basin)
    )
    lake_mask = np.zeros(grid.n, dtype=bool)
    if lake_candidates.any():
        count = min(max(4, int(0.00006 * grid.n)),
                    int(np.count_nonzero(lake_candidates)))
        score = accumulation + np.where(basin, 80.0, 0.0) - np.maximum(rel, 0.0) / 9.0
        chosen = np.flatnonzero(lake_candidates)[
            np.argsort(score[lake_candidates])[-count:]
        ]
        lake_mask[chosen] = True
    if lake_mask.sum() < max(2, int(0.000035 * grid.n)):
        basin_candidates = finite_land & basin & (distance >= 4.0) & (rel < 700.0)
        if basin_candidates.any():
            score = accumulation + (700.0 - np.clip(rel, 0.0, 700.0)) / 12.0
            count = min(max(3, int(0.00005 * grid.n)),
                        int(np.count_nonzero(basin_candidates)))
            chosen = np.flatnonzero(basin_candidates)[
                np.argsort(score[basin_candidates])[-count:]
            ]
            lake_mask[chosen] = True
    if lake_mask.any():
        lake_domain = finite_land & (distance >= 3.0) & ((rel < 760.0) | basin)
        lake_mask = (_dilate_mask(grid, lake_mask, passes=1) & lake_domain) | lake_mask
    lake_basin_rank = np.zeros(grid.n, dtype=np.float64)
    if lake_mask.any():
        lake_domain = finite_land & (distance >= 2.0) & ((rel < 900.0) | basin)
        lake_basin_rank = _ranked_halo(
            grid,
            lake_mask,
            domain=lake_domain,
            passes=2,
            decay=0.52,
        )

    meander_scroll_rank = np.clip(
        meander_belt_rank
        * np.maximum(floodplain_rank, 0.36)
        * finite_land
        * (rel < 980.0)
        * (distance >= 1.0),
        0.0,
        1.0,
    )
    meander_scroll_rank[meander_scroll_rank < 0.10] = 0.0
    floodplain_swale_rank = np.clip(
        floodplain_rank
        * (0.55 + 0.45 * basin.astype(np.float64))
        * finite_land
        * (rel < 820.0)
        * (river_path_rank < 0.72),
        0.0,
        1.0,
    )
    floodplain_swale_rank[floodplain_swale_rank < 0.10] = 0.0
    lake_shoreline_rank = np.zeros(grid.n, dtype=np.float64)
    if lake_mask.any():
        lake_shore_domain = finite_land & (lake_basin_rank > 0.0) & ~lake_mask
        lake_rim_seed = _dilate_mask(grid, lake_mask, passes=1) & lake_shore_domain
        if lake_rim_seed.any():
            lake_shoreline_rank = _ranked_halo(
                grid,
                lake_rim_seed,
                domain=lake_shore_domain,
                passes=1,
                decay=0.56,
            )
            lake_shoreline_rank *= (
                0.56
                + 0.28 * np.clip(lake_basin_rank, 0.0, 1.0)
                + 0.16 * np.clip((900.0 - rel) / 900.0, 0.0, 1.0)
            )
        lake_shoreline_rank[~lake_shore_domain] = 0.0
        lake_shoreline_rank[lake_shoreline_rank < 0.12] = 0.0

    river_carve_rank = np.maximum(river_path_rank, river_rank * 0.38)
    coastal_acc_cut = (
        float(np.nanpercentile(accumulation[coastal_land], 88.0))
        if np.any(coastal_land)
        else np.inf
    )
    delta_mouth = (
        coastal_land
        & (river_path_rank >= 0.62)
        & (accumulation >= coastal_acc_cut)
    )
    ocean_coast_distance = _ocean_coast_distance(grid, land, ocean, max_passes=4)
    delta_seed = _ocean_neighbors_of(grid, delta_mouth, ocean)
    delta_domain = ocean & (ocean_coast_distance >= 0) & (ocean_coast_distance <= 3) & (
        rel > -1800.0
    )
    delta_fan_rank = _ranked_halo(
        grid,
        delta_seed,
        domain=delta_domain,
        passes=2,
        decay=0.62,
    )
    delta_mask = delta_fan_rank > 0.18
    delta_plain_domain = finite_land & (distance <= 2.0) & (rel < 520.0)
    delta_plain_mask = _ranked_halo(
        grid,
        delta_mouth,
        domain=delta_plain_domain,
        passes=2,
        decay=0.55,
    ) > 0.18

    hydrology_delta = np.zeros(grid.n, dtype=np.float64)
    river_cells = finite_land & (river_carve_rank > 0.08)
    hydrology_delta[river_cells] -= 18.0 + 125.0 * river_carve_rank[river_cells]
    floodplain_cells = finite_land & (floodplain_rank > 0.0) & (river_carve_rank <= 0.18)
    hydrology_delta[floodplain_cells] -= 8.0 + 28.0 * floodplain_rank[floodplain_cells]
    lake_basin_cells = finite_land & (lake_basin_rank > 0.0)
    hydrology_delta[lake_basin_cells] -= np.minimum(
        0.30 * np.maximum(rel[lake_basin_cells], 0.0),
        24.0 + 82.0 * lake_basin_rank[lake_basin_cells],
    )
    hydrology_delta[lake_mask] -= np.minimum(0.40 * np.maximum(rel[lake_mask], 0.0),
                                             150.0)
    hydrology_delta[delta_plain_mask] -= 10.0 + 34.0 * np.maximum(
        river_path_rank[delta_plain_mask],
        floodplain_rank[delta_plain_mask],
    )
    hydrology_delta[delta_mouth] -= 18.0 + 45.0 * river_path_rank[delta_mouth]
    hydrology_delta[delta_mask] += 38.0 + 82.0 * delta_fan_rank[delta_mask]

    base_hydrology_delta = hydrology_delta.copy()
    lowland_gate = np.clip((1120.0 - rel) / 1120.0, 0.0, 1.0) * finite_land
    scroll_texture = _signed_texture(grid, context.detail_seed + 1703, scale=5)
    swale_texture = _ridged_texture(grid, context.detail_seed + 1711)
    lake_texture = _signed_texture(grid, context.detail_seed + 1721, scale=4)
    levee_seed = finite_land & (river_path_rank >= 0.30)
    levee_domain = (
        finite_land
        & (floodplain_rank > 0.0)
        & (river_path_rank < 0.30)
        & (rel < 980.0)
    )
    natural_levee_rank = _ranked_halo(
        grid,
        levee_seed,
        domain=levee_domain,
        passes=1,
        decay=0.55,
    )
    scroll_support = np.sqrt(np.clip(meander_scroll_rank, 0.0, 1.0))
    swale_support = np.sqrt(np.clip(floodplain_swale_rank, 0.0, 1.0))
    levee_support = np.sqrt(np.clip(natural_levee_rank, 0.0, 1.0))
    lake_shore_support = np.sqrt(np.clip(lake_shoreline_rank, 0.0, 1.0))
    floodplain_support = np.sqrt(
        np.clip(np.maximum(floodplain_rank, meander_belt_rank), 0.0, 1.0)
    )
    scroll_delta = (
        scroll_support
        * lowland_gate
        * (13.0 * scroll_texture + 8.0 * (swale_texture - 0.42))
    )
    swale_delta = (
        -swale_support
        * lowland_gate
        * (4.0 + 13.0 * (0.65 + 0.35 * swale_texture))
    )
    levee_delta = (
        levee_support
        * lowland_gate
        * (6.0 + 17.0 * floodplain_support)
    )
    lake_shore_delta = (
        lake_shore_support
        * lowland_gate
        * (3.5 + 7.5 * np.clip(0.55 + 0.45 * lake_texture, 0.0, 1.0))
    )
    fluvial_microrelief_delta = np.clip(
        scroll_delta + swale_delta + levee_delta + lake_shore_delta,
        -40.0,
        34.0,
    )
    fluvial_microrelief_delta *= _polar_detail_damping(grid)
    fluvial_microrelief_delta[~finite_land] = 0.0
    hydrology_delta += fluvial_microrelief_delta

    baseline_refined = rel + base_hydrology_delta
    baseline_refined[land] = np.maximum(baseline_refined[land], 2.0)
    baseline_refined[ocean] = np.minimum(baseline_refined[ocean], -2.0)
    refined = rel + hydrology_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean] = np.minimum(refined[ocean], -2.0)
    applied_microrelief_delta = refined - baseline_refined
    lowland_baseline = refined.copy()
    (
        lowland_alluvial_delta,
        alluvial_fan_rank,
        lowland_plain_rank,
        piedmont_apron_rank,
    ) = _lowland_alluvial_microrelief_delta(
        grid,
        rel,
        refined,
        finite_land,
        distance,
        detail_region,
        hierarchy,
        shoulder,
        apron,
        basin,
        river_path_rank,
        basin_trunk_rank,
        floodplain_rank,
        lake_basin_rank,
        delta_plain_mask,
        context.detail_seed,
    )
    hydrology_delta += lowland_alluvial_delta
    refined = rel + hydrology_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean] = np.minimum(refined[ocean], -2.0)
    lowland_alluvial_delta = refined - lowland_baseline

    distance_out = np.where(finite_land, distance, -1.0)
    return _Microgeomorphology(
        refined_rel=refined,
        hydrology_delta=refined - rel,
        fluvial_microrelief_delta=applied_microrelief_delta,
        lowland_alluvial_microrelief_delta=lowland_alluvial_delta,
        flow_accumulation=accumulation,
        river_rank=river_rank,
        river_path_rank=river_path_rank,
        basin_trunk_rank=basin_trunk_rank,
        floodplain_rank=floodplain_rank,
        meander_belt_rank=meander_belt_rank,
        meander_scroll_rank=meander_scroll_rank,
        floodplain_swale_rank=floodplain_swale_rank,
        alluvial_fan_rank=alluvial_fan_rank,
        lowland_plain_rank=lowland_plain_rank,
        piedmont_apron_rank=piedmont_apron_rank,
        lake_basin_rank=lake_basin_rank,
        lake_shoreline_rank=lake_shoreline_rank,
        delta_fan_rank=delta_fan_rank,
        drainage_basin_id=basin_id,
        river_receiver=receiver,
        lake_mask=lake_mask,
        delta_mask=delta_mask,
        delta_plain_mask=delta_plain_mask,
        land_coast_distance=distance_out,
    )


def _lowland_alluvial_microrelief_delta(
    grid: SphereGrid,
    rel: np.ndarray,
    refined_rel: np.ndarray,
    finite_land: np.ndarray,
    land_distance: np.ndarray,
    detail_region: np.ndarray,
    hierarchy: np.ndarray,
    shoulder: np.ndarray,
    apron: np.ndarray,
    basin: np.ndarray,
    river_path_rank: np.ndarray,
    basin_trunk_rank: np.ndarray,
    floodplain_rank: np.ndarray,
    lake_basin_rank: np.ndarray,
    delta_plain_mask: np.ndarray,
    detail_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Add subtle lowland plain, fan, and piedmont relief on land only."""
    rel = np.asarray(rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    finite_land = np.asarray(finite_land, dtype=bool)
    distance = np.asarray(land_distance, dtype=np.float64)
    detail_region = np.asarray(detail_region)
    hierarchy = np.asarray(hierarchy)
    shoulder = np.asarray(shoulder)
    apron = np.asarray(apron)
    basin = np.asarray(basin, dtype=bool)
    river_path_rank = np.asarray(river_path_rank, dtype=np.float64)
    basin_trunk_rank = np.asarray(basin_trunk_rank, dtype=np.float64)
    floodplain_rank = np.asarray(floodplain_rank, dtype=np.float64)
    lake_basin_rank = np.asarray(lake_basin_rank, dtype=np.float64)
    delta_plain_mask = np.asarray(delta_plain_mask, dtype=bool)

    delta = np.zeros(rel.shape, dtype=np.float64)
    fan_rank = np.zeros(rel.shape, dtype=np.float64)
    plain_rank = np.zeros(rel.shape, dtype=np.float64)
    piedmont_rank = np.zeros(rel.shape, dtype=np.float64)
    if not np.any(finite_land):
        return delta, fan_rank, plain_rank, piedmont_rank

    low_polar = np.abs(grid.lat) < 62.0
    lowland_gate = np.clip((950.0 - refined_rel) / 950.0, 0.0, 1.0)
    midland_gate = np.clip((1700.0 - refined_rel) / 1700.0, 0.0, 1.0)
    hydro_corridor = np.maximum.reduce([
        river_path_rank,
        basin_trunk_rank,
        floodplain_rank,
        0.55 * lake_basin_rank,
    ])
    mountain_seed = (
        finite_land
        & low_polar
        & ((hierarchy > 0) | (shoulder > 0) | (apron > 0))
        & (refined_rel > 420.0)
    )
    mountain_front = _ranked_halo(
        grid,
        mountain_seed,
        domain=finite_land & low_polar & (refined_rel < 1800.0) & (distance >= 1.0),
        passes=3,
        decay=0.54,
    )
    piedmont_domain = (
        finite_land
        & low_polar
        & (distance >= 2.0)
        & (refined_rel > 120.0)
        & (refined_rel < 1550.0)
        & (mountain_front > 0.18)
        & (hierarchy <= 1)
    )
    fan_seed = (
        piedmont_domain
        & (hydro_corridor > 0.12)
        & (floodplain_rank < 0.78)
        & (lake_basin_rank < 0.72)
    )
    if np.any(fan_seed):
        fan_rank = _ranked_halo(
            grid,
            fan_seed,
            domain=piedmont_domain,
            passes=2,
            decay=0.48,
        )
        fan_rank *= np.clip(
            0.38
            + 0.38 * mountain_front
            + 0.24 * np.sqrt(np.clip(hydro_corridor, 0.0, 1.0)),
            0.0,
            1.0,
        )
        fan_rank[fan_rank < 0.10] = 0.0

    plain_domain = (
        finite_land
        & low_polar
        & (distance >= 1.0)
        & (refined_rel < 720.0)
        & (hierarchy <= 0)
        & (mountain_front < 0.58)
        & (
            basin
            | delta_plain_mask
            | (floodplain_rank > 0.06)
            | (lake_basin_rank > 0.10)
            | (hydro_corridor > 0.18)
        )
    )
    if np.any(plain_domain):
        plain_seed = (
            plain_domain
            & (
                basin
                | delta_plain_mask
                | (floodplain_rank > 0.10)
                | (lake_basin_rank > 0.12)
                | ((hydro_corridor > 0.20) & (refined_rel < 620.0))
            )
        )
        if np.any(plain_seed):
            plain_rank = _ranked_halo(
                grid,
                plain_seed,
                domain=plain_domain,
                passes=2,
                decay=0.50,
            )
        plain_rank = np.maximum(
            plain_rank,
            np.where(
                plain_domain,
                np.clip(
                    0.34 * basin.astype(np.float64)
                    + 0.24 * floodplain_rank
                    + 0.16 * lake_basin_rank
                    + 0.20 * delta_plain_mask.astype(np.float64)
                    + 0.12 * np.sqrt(np.clip(hydro_corridor, 0.0, 1.0)),
                    0.0,
                    0.78,
                ),
                0.0,
            ),
        )
        plain_rank[plain_rank < 0.12] = 0.0

    if np.any(piedmont_domain):
        piedmont_rank = np.clip(
            mountain_front
            * midland_gate
            * (0.42 + 0.32 * fan_rank + 0.26 * np.sqrt(np.clip(hydro_corridor, 0.0, 1.0))),
            0.0,
            1.0,
        )
        piedmont_rank[piedmont_rank < 0.12] = 0.0
        piedmont_rank = np.maximum(piedmont_rank, 0.55 * fan_rank)

    local_mean = _local_domain_neighbor_mean(grid, refined_rel, finite_land)
    texture = _signed_texture(grid, int(detail_seed) + 5303, scale=5)
    ridged = _ridged_texture(grid, int(detail_seed) + 5323)
    plain_flatten = np.clip(local_mean - refined_rel, -26.0, 26.0)
    plain_delta = (
        plain_rank
        * lowland_gate
        * (0.40 * plain_flatten + 2.2 * texture)
    )
    fan_delta = (
        fan_rank
        * midland_gate
        * (5.0 + 18.0 * ridged - 7.0 * np.sqrt(np.clip(river_path_rank, 0.0, 1.0)))
    )
    piedmont_delta = (
        piedmont_rank
        * midland_gate
        * (0.22 * np.clip(local_mean - refined_rel, -32.0, 32.0) + 3.4 * texture)
    )
    delta = np.clip(plain_delta + fan_delta + piedmont_delta, -28.0, 34.0)
    delta *= _polar_detail_damping(grid)
    delta[~finite_land] = 0.0
    return (
        delta,
        np.clip(fan_rank, 0.0, 1.0),
        np.clip(plain_rank, 0.0, 1.0),
        np.clip(piedmont_rank, 0.0, 1.0),
    )


@dataclass
class _CoastalMorphology:
    refined_rel: np.ndarray
    coastal_delta: np.ndarray
    coastal_process_microrelief_delta: np.ndarray
    coastal_depositional_microrelief_delta: np.ndarray
    coastal_plain_rank: np.ndarray
    coastal_cliff_rank: np.ndarray
    shoreface_rank: np.ndarray
    barrier_lagoon_rank: np.ndarray
    estuary_rank: np.ndarray
    delta_distributary_rank: np.ndarray
    estuary_funnel_rank: np.ndarray
    barrier_spit_rank: np.ndarray
    delta_mouth_bar_rank: np.ndarray
    estuary_tidal_channel_rank: np.ndarray
    coastal_depositional_plain_rank: np.ndarray
    strandplain_rank: np.ndarray
    tidal_flat_rank: np.ndarray


def _selected_snapshot_coastal_morphology(
    context: _RefinementContext,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    micro: _Microgeomorphology,
) -> _CoastalMorphology:
    """Add topology-preserving coastal plain, cliff, shoreface, and estuary ranks."""
    grid = context.target_grid
    rel = np.asarray(refined_rel, dtype=np.float64)
    land = rel >= 0.0
    ocean = ~land
    zeros = np.zeros(grid.n, dtype=np.float64)
    if not land.any() or not ocean.any():
        return _CoastalMorphology(
            refined_rel=rel.copy(),
            coastal_delta=zeros,
            coastal_process_microrelief_delta=zeros,
            coastal_depositional_microrelief_delta=zeros,
            coastal_plain_rank=zeros,
            coastal_cliff_rank=zeros,
            shoreface_rank=zeros,
            barrier_lagoon_rank=zeros,
            estuary_rank=zeros,
            delta_distributary_rank=zeros,
            estuary_funnel_rank=zeros,
            barrier_spit_rank=zeros,
            delta_mouth_bar_rank=zeros,
            estuary_tidal_channel_rank=zeros,
            coastal_depositional_plain_rank=zeros,
            strandplain_rank=zeros,
            tidal_flat_rank=zeros,
        )

    land_distance = np.asarray(micro.land_coast_distance, dtype=np.float64)
    ocean_distance = np.asarray(marine.ocean_coast_distance, dtype=np.float64)
    detail_region = _int_field(
        context, "terrain.continental_detail_region_code", default=0)
    hierarchy = _int_field(context, "terrain.orogenic_parent_hierarchy", default=0)
    active_margin = _mask(context, "boundary__active_margin") | _mask(
        context, "object__boundary_active_margin")
    passive_margin = _mask(context, "boundary__passive_margin") | _mask(
        context, "object__boundary_passive_margin")
    trench = _mask(context, "boundary__trench") | _mask(
        context, "object__margin_trench")
    ridge_texture = _ridged_texture(grid, context.detail_seed + 2309)
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    low_polar = lat_abs <= 72.0

    platform = detail_region == CONT_DETAIL_PLATFORM
    basin = np.isin(detail_region, [CONT_DETAIL_BASIN, CONT_DETAIL_RIFT_BASIN])
    orogen = (detail_region == CONT_DETAIL_OROGEN) | (hierarchy > 0)
    land_coast = land & (land_distance >= 0.0)
    land_falloff = np.exp(-np.maximum(land_distance, 0.0) / 2.15)
    low_relief = np.clip((620.0 - rel) / 620.0, 0.0, 1.0)
    plain_support = (
        0.32
        + 0.24 * passive_margin.astype(np.float64)
        + 0.20 * platform.astype(np.float64)
        + 0.18 * basin.astype(np.float64)
    )
    coastal_plain_rank = (
        plain_support
        * land_falloff
        * low_relief
        * land_coast
        * (land_distance <= 2.5)
        * (rel < 560.0)
        * low_polar
        * (passive_margin | platform | basin)
        * ~active_margin
        * ~trench
        * ~orogen
    )

    cliff_support = (
        0.36
        + 0.28 * active_margin.astype(np.float64)
        + 0.22 * trench.astype(np.float64)
        + 0.22 * orogen.astype(np.float64)
    )
    coastal_cliff_rank = (
        cliff_support
        * np.exp(-np.maximum(land_distance, 0.0) / 0.78)
        * np.clip((rel - 95.0) / 980.0, 0.0, 1.0)
        * land_coast
        * (land_distance <= 0.5)
        * (rel > 180.0)
        * low_polar
        * (active_margin | trench | orogen)
    )
    coastal_cliff_rank = np.clip(coastal_cliff_rank, 0.0, 1.0)
    coastal_plain_rank = np.clip(coastal_plain_rank, 0.0, 1.0)
    coastal_plain_rank[coastal_plain_rank < 0.055] = 0.0
    coastal_cliff_rank[coastal_cliff_rank < 0.16] = 0.0

    ocean_coast = ocean & (ocean_distance >= 0.0)
    shoreface_domain = (
        ocean_coast
        & (ocean_distance <= 1.65)
        & (rel < -3.0)
        & (rel > -720.0)
        & low_polar
        & ~trench
    )
    depth_factor = np.clip((1050.0 + rel) / 1050.0, 0.0, 1.0)
    shoreface_rank = (
        shoreface_domain
        * np.clip((2.75 - ocean_distance) / 2.75, 0.0, 1.0)
        * depth_factor
        * (0.52 + 0.22 * passive_margin.astype(np.float64)
           + 0.10 * active_margin.astype(np.float64))
    )
    shoreface_rank = np.clip(shoreface_rank, 0.0, 1.0)
    shoreface_rank[shoreface_rank < 0.12] = 0.0

    passive_like = passive_margin | (~active_margin & ~trench & (shoreface_rank > 0.16))
    barrier_domain = (
        shoreface_domain
        & passive_like
        & ~active_margin
        & (lat_abs <= 58.0)
        & (ocean_distance <= 1.25)
        & (rel > -290.0)
        & (rel < -8.0)
    )
    barrier_score = (
        0.44 * shoreface_rank
        + 0.22 * np.clip((360.0 + rel) / 360.0, 0.0, 1.0)
        + 0.20 * ridge_texture
        + 0.14 * np.clip((58.0 - lat_abs) / 58.0, 0.0, 1.0)
    )
    barrier_lagoon_rank = np.zeros(grid.n, dtype=np.float64)
    if barrier_domain.any():
        cut = max(0.60, float(np.nanpercentile(barrier_score[barrier_domain], 92.0)))
        barrier_seed = barrier_domain & (barrier_score >= cut)
        barrier_lagoon_rank = _ranked_halo(
            grid,
            barrier_seed,
            domain=barrier_domain,
            passes=1,
            decay=0.46,
        )
        barrier_lagoon_rank *= np.maximum(barrier_score, 0.58)
        barrier_lagoon_rank[~barrier_domain] = 0.0
        barrier_lagoon_rank[barrier_lagoon_rank < 0.18] = 0.0

    land_neighbor_fraction = _neighbor_mask_fraction(grid, land)
    coastal_confinement = np.clip((land_neighbor_fraction - 0.16) / 0.42, 0.0, 1.0)
    open_shelf_context = np.clip((0.48 - land_neighbor_fraction) / 0.36, 0.0, 1.0)
    river_mouth_rank = np.maximum(
        np.asarray(micro.delta_fan_rank, dtype=np.float64),
        np.asarray(micro.delta_mask, dtype=np.float64) * 0.65,
    )
    shallow_delta_gate = np.clip((760.0 + rel) / 760.0, 0.0, 1.0)
    estuary_depth_gate = np.clip((1120.0 + rel) / 1120.0, 0.0, 1.0)
    delta_context = np.clip(
        0.42
        + 0.46 * open_shelf_context
        + 0.12 * passive_margin.astype(np.float64)
        - 0.18 * active_margin.astype(np.float64)
        - 0.24 * coastal_confinement,
        0.0,
        1.0,
    )
    estuary_context = np.clip(
        0.20
        + 0.64 * coastal_confinement
        + 0.12 * active_margin.astype(np.float64)
        + 0.10 * shoreface_rank
        - 0.20 * open_shelf_context
        - 0.18 * barrier_lagoon_rank,
        0.0,
        1.0,
    )
    estuary_domain = (
        ocean_coast
        & (ocean_distance <= 3.0)
        & (rel > -950.0)
        & low_polar
        & ~trench
    )
    estuary_rank = (
        river_mouth_rank
        * estuary_context
        * estuary_depth_gate
        * (1.0 - 0.42 * open_shelf_context)
    )
    estuary_rank = np.where(estuary_domain, estuary_rank, 0.0)
    estuary_rank = np.clip(estuary_rank, 0.0, 1.0)

    delta_distributary_rank = np.clip(
        river_mouth_rank
        * delta_context
        * shallow_delta_gate
        * (ocean_distance <= 3.0)
        * ocean
        * low_polar
        * ~trench,
        0.0,
        1.0,
    )
    delta_distributary_rank[delta_distributary_rank < 0.16] = 0.0
    estuary_funnel_rank = np.clip(
        estuary_rank
        * np.clip((3.35 - ocean_distance) / 3.35, 0.0, 1.0)
        * (rel > -1050.0)
        * ocean
        * low_polar,
        0.0,
        1.0,
    )
    estuary_funnel_rank[estuary_funnel_rank < 0.12] = 0.0
    barrier_spit_rank = np.clip(
        barrier_lagoon_rank
        * np.clip((1.55 - ocean_distance) / 1.55, 0.0, 1.0)
        * ocean
        * (lat_abs <= 58.0)
        * ~trench,
        0.0,
        1.0,
    )
    barrier_spit_rank[barrier_spit_rank < 0.16] = 0.0
    delta_mouth_bar_rank = np.zeros(grid.n, dtype=np.float64)
    delta_bar_score = (
        delta_distributary_rank
        * open_shelf_context
        * shallow_delta_gate
        * (1.0 - 0.42 * np.clip(estuary_funnel_rank, 0.0, 1.0))
        * (1.0 - 0.36 * np.clip(barrier_spit_rank, 0.0, 1.0))
    )
    delta_bar_domain = (
        ocean
        & (delta_distributary_rank > 0.0)
        & (ocean_distance >= 0.0)
        & (ocean_distance <= 2.0)
        & (rel > -820.0)
        & (delta_distributary_rank >= 0.82 * estuary_funnel_rank)
        & (barrier_spit_rank < 0.46)
    )
    if np.any(delta_bar_domain):
        cut = max(0.13, float(np.nanpercentile(delta_bar_score[delta_bar_domain], 68.0)))
        delta_bar_seed = delta_bar_domain & (delta_bar_score >= cut)
        if np.any(delta_bar_seed):
            delta_mouth_bar_rank = _ranked_halo(
                grid,
                delta_bar_seed,
                domain=delta_bar_domain,
                passes=1,
                decay=0.34,
            )
            delta_mouth_bar_rank *= np.clip(delta_bar_score / max(cut, 1.0e-9), 0.30, 1.0)
            delta_mouth_bar_rank[~delta_bar_domain] = 0.0
            delta_mouth_bar_rank[delta_mouth_bar_rank < 0.12] = 0.0

    estuary_tidal_channel_rank = np.zeros(grid.n, dtype=np.float64)
    tidal_channel_score = (
        estuary_funnel_rank
        * coastal_confinement
        * estuary_depth_gate
        * (1.0 - 0.52 * np.clip(delta_distributary_rank, 0.0, 1.0))
        * (1.0 - 0.38 * np.clip(barrier_spit_rank, 0.0, 1.0))
    )
    tidal_channel_domain = (
        ocean
        & (estuary_funnel_rank > 0.0)
        & (ocean_distance >= 0.0)
        & (ocean_distance <= 2.8)
        & (rel > -1050.0)
        & (estuary_funnel_rank >= 0.72 * delta_distributary_rank)
        & (barrier_spit_rank < 0.42)
    )
    if np.any(tidal_channel_domain):
        cut = max(0.10, float(np.nanpercentile(tidal_channel_score[tidal_channel_domain], 62.0)))
        tidal_channel_seed = tidal_channel_domain & (tidal_channel_score >= cut)
        if np.any(tidal_channel_seed):
            estuary_tidal_channel_rank = _ranked_halo(
                grid,
                tidal_channel_seed,
                domain=tidal_channel_domain,
                passes=1,
                decay=0.24,
            )
            estuary_tidal_channel_rank *= np.clip(
                tidal_channel_score / max(cut, 1.0e-9), 0.28, 1.0)
            estuary_tidal_channel_rank[~tidal_channel_domain] = 0.0
            estuary_tidal_channel_rank[estuary_tidal_channel_rank < 0.10] = 0.0

    coastal_delta = np.zeros(grid.n, dtype=np.float64)
    plain_cells = land & (coastal_plain_rank > 0.05)
    plain_target = 16.0 + 38.0 * np.maximum(land_distance, 0.0)
    coastal_delta[plain_cells] -= np.minimum(
        46.0 * coastal_plain_rank[plain_cells],
        np.maximum(0.0, rel[plain_cells] - plain_target[plain_cells]) * 0.18,
    )

    cliff_cells = land & (coastal_cliff_rank > 0.08)
    coastal_delta[cliff_cells] += 18.0 * coastal_cliff_rank[cliff_cells]

    d = np.maximum(ocean_distance, 0.0)
    shoreface_target = -(28.0 + 68.0 * d + 34.0 * d ** 1.28)
    shore_cells = ocean & (shoreface_rank > 0.05)
    coastal_delta[shore_cells] += np.clip(
        0.20 * (shoreface_target[shore_cells] - rel[shore_cells]),
        -34.0,
        118.0,
    ) * shoreface_rank[shore_cells]

    barrier_cells = ocean & (barrier_lagoon_rank > 0.05)
    barrier_target = -13.0 - 13.0 * (1.0 - barrier_lagoon_rank)
    coastal_delta[barrier_cells] += np.minimum(
        62.0 * barrier_lagoon_rank[barrier_cells],
        np.maximum(0.0, barrier_target[barrier_cells] - rel[barrier_cells]) * 0.48,
    )

    estuary_cells = ocean & (estuary_rank > 0.05)
    estuary_target = -42.0 - 48.0 * estuary_rank
    coastal_delta[estuary_cells] -= np.minimum(
        42.0 * estuary_rank[estuary_cells],
        np.maximum(0.0, rel[estuary_cells] - estuary_target[estuary_cells]) * 0.35,
    )

    cliff_ocean_seed = _ocean_neighbors_of(grid, coastal_cliff_rank > 0.18, ocean)
    if cliff_ocean_seed.any():
        cliff_ocean_rank = _ranked_halo(
            grid,
            cliff_ocean_seed,
            domain=ocean & (ocean_distance >= 0.0) & (ocean_distance <= 1.5) & ~trench,
            passes=1,
            decay=0.52,
        )
        coastal_delta[ocean & (cliff_ocean_rank > 0.0)] -= 30.0 * cliff_ocean_rank[
            ocean & (cliff_ocean_rank > 0.0)
        ]

    base_coastal_delta = coastal_delta.copy()
    process_texture = _signed_texture(grid, context.detail_seed + 2503, scale=5)
    process_ridge = _ridged_texture(grid, context.detail_seed + 2511)
    shallow_gate = np.clip((1180.0 + rel) / 1180.0, 0.0, 1.0) * ocean
    delta_support = np.sqrt(np.clip(delta_distributary_rank, 0.0, 1.0))
    estuary_support = np.sqrt(np.clip(estuary_funnel_rank, 0.0, 1.0))
    barrier_support = np.sqrt(np.clip(barrier_spit_rank, 0.0, 1.0))
    delta_process_delta = (
        delta_support
        * shallow_gate
        * (7.0 + 22.0 * process_ridge + 6.0 * np.maximum(process_texture, 0.0))
    )
    estuary_process_delta = (
        -estuary_support
        * (1.0 - 0.86 * np.clip(delta_support, 0.0, 1.0))
        * shallow_gate
        * (8.0 + 26.0 * (0.55 + 0.45 * process_ridge))
    )
    barrier_process_delta = (
        barrier_support
        * shallow_gate
        * (3.0 + 9.0 * (0.58 + 0.42 * process_ridge))
    )
    mouth_bar_delta = (
        delta_mouth_bar_rank
        * shallow_gate
        * (2.5 + 7.5 * (0.50 + 0.50 * process_ridge))
    )
    tidal_channel_delta = (
        -estuary_tidal_channel_rank
        * shallow_gate
        * (3.5 + 10.0 * (0.45 + 0.55 * process_ridge))
    )
    coastal_process_microrelief_delta = np.clip(
        delta_process_delta
        + estuary_process_delta
        + barrier_process_delta
        + mouth_bar_delta
        + tidal_channel_delta,
        -34.0,
        30.0,
    )
    delta_mouth_floor = (
        delta_support
        * shallow_gate
        * (3.0 + 8.0 * process_ridge)
    )
    delta_mouth = delta_support > 0.35
    coastal_process_microrelief_delta[delta_mouth] = np.maximum(
        coastal_process_microrelief_delta[delta_mouth],
        delta_mouth_floor[delta_mouth],
    )
    coastal_process_microrelief_delta *= _polar_detail_damping(grid)
    coastal_process_microrelief_delta[~ocean] = 0.0
    coastal_delta += coastal_process_microrelief_delta

    process_baseline = rel + base_coastal_delta
    process_baseline[land] = np.maximum(process_baseline[land], 2.0)
    process_baseline[ocean] = np.minimum(process_baseline[ocean], -2.0)
    process_refined = rel + coastal_delta
    process_refined[land] = np.maximum(process_refined[land], 2.0)
    process_refined[ocean] = np.minimum(process_refined[ocean], -2.0)
    applied_process_delta = process_refined - process_baseline
    (
        coastal_depositional_delta,
        coastal_depositional_plain_rank,
        strandplain_rank,
        tidal_flat_rank,
    ) = _coastal_depositional_microrelief_delta(
        grid,
        rel,
        process_refined,
        land,
        land_distance,
        active_margin,
        trench,
        orogen,
        coastal_plain_rank,
        shoreface_rank,
        barrier_lagoon_rank,
        estuary_rank,
        np.asarray(micro.delta_plain_mask, dtype=bool),
        np.asarray(micro.lowland_plain_rank, dtype=np.float64),
        np.asarray(micro.floodplain_rank, dtype=np.float64),
        context.detail_seed,
    )
    coastal_delta += coastal_depositional_delta
    refined = rel + coastal_delta
    refined[land] = np.maximum(refined[land], 2.0)
    refined[ocean] = np.minimum(refined[ocean], -2.0)
    coastal_depositional_delta = refined - process_refined
    return _CoastalMorphology(
        refined_rel=refined,
        coastal_delta=refined - rel,
        coastal_process_microrelief_delta=applied_process_delta,
        coastal_depositional_microrelief_delta=coastal_depositional_delta,
        coastal_plain_rank=coastal_plain_rank,
        coastal_cliff_rank=coastal_cliff_rank,
        shoreface_rank=shoreface_rank,
        barrier_lagoon_rank=np.clip(barrier_lagoon_rank, 0.0, 1.0),
        estuary_rank=estuary_rank,
        delta_distributary_rank=delta_distributary_rank,
        estuary_funnel_rank=estuary_funnel_rank,
        barrier_spit_rank=barrier_spit_rank,
        delta_mouth_bar_rank=np.clip(delta_mouth_bar_rank, 0.0, 1.0),
        estuary_tidal_channel_rank=np.clip(estuary_tidal_channel_rank, 0.0, 1.0),
        coastal_depositional_plain_rank=coastal_depositional_plain_rank,
        strandplain_rank=strandplain_rank,
        tidal_flat_rank=tidal_flat_rank,
    )


def _coastal_depositional_microrelief_delta(
    grid: SphereGrid,
    rel: np.ndarray,
    refined_rel: np.ndarray,
    land: np.ndarray,
    land_distance: np.ndarray,
    active_margin: np.ndarray,
    trench: np.ndarray,
    orogen: np.ndarray,
    coastal_plain_rank: np.ndarray,
    shoreface_rank: np.ndarray,
    barrier_lagoon_rank: np.ndarray,
    estuary_rank: np.ndarray,
    delta_plain_mask: np.ndarray,
    lowland_plain_rank: np.ndarray,
    floodplain_rank: np.ndarray,
    detail_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Add land-side coastal depositional microrelief without changing coastline."""
    rel = np.asarray(rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    land = np.asarray(land, dtype=bool)
    distance = np.asarray(land_distance, dtype=np.float64)
    active_margin = np.asarray(active_margin, dtype=bool)
    trench = np.asarray(trench, dtype=bool)
    orogen = np.asarray(orogen, dtype=bool)
    coastal_plain_rank = np.asarray(coastal_plain_rank, dtype=np.float64)
    shoreface_rank = np.asarray(shoreface_rank, dtype=np.float64)
    barrier_lagoon_rank = np.asarray(barrier_lagoon_rank, dtype=np.float64)
    estuary_rank = np.asarray(estuary_rank, dtype=np.float64)
    delta_plain_mask = np.asarray(delta_plain_mask, dtype=bool)
    lowland_plain_rank = np.asarray(lowland_plain_rank, dtype=np.float64)
    floodplain_rank = np.asarray(floodplain_rank, dtype=np.float64)
    ocean = ~land
    low_polar = np.abs(grid.lat) <= 58.0
    low_land = (
        land
        & low_polar
        & (distance >= 0.0)
        & (distance <= 2.6)
        & (refined_rel >= 2.0)
        & (refined_rel < 320.0)
        & ~active_margin
        & ~trench
        & ~orogen
    )
    delta = np.zeros(rel.shape, dtype=np.float64)
    depositional_rank = np.zeros(rel.shape, dtype=np.float64)
    strandplain_rank = np.zeros(rel.shape, dtype=np.float64)
    tidal_flat_rank = np.zeros(rel.shape, dtype=np.float64)
    if not np.any(low_land):
        return delta, depositional_rank, strandplain_rank, tidal_flat_rank

    shore_seed = _land_neighbors_of(
        grid,
        ocean & (shoreface_rank > 0.18),
        land,
    )
    barrier_seed = _land_neighbors_of(
        grid,
        ocean & (barrier_lagoon_rank > 0.16),
        land,
    )
    estuary_seed = _land_neighbors_of(
        grid,
        ocean & (estuary_rank > 0.12),
        land,
    )
    plain_lowland_seed = (
        (coastal_plain_rank > 0.055)
        & (distance <= 2.0)
        & (
            (lowland_plain_rank > 0.12)
            | (floodplain_rank > 0.10)
            | (distance <= 1.2)
        )
    )
    shore_lowland_seed = (
        shore_seed
        & (distance <= 1.35)
        & (
            (lowland_plain_rank > 0.20)
            | (floodplain_rank > 0.16)
            | (coastal_plain_rank > 0.07)
        )
    )
    process_seed = (
        barrier_seed
        | estuary_seed
        | delta_plain_mask
        | plain_lowland_seed
        | shore_lowland_seed
    ) & low_land
    if not np.any(process_seed):
        return delta, depositional_rank, strandplain_rank, tidal_flat_rank

    process_halo = _dilate_mask(grid, process_seed, passes=1)
    supported_lowland = (
        (distance <= 1.8)
        & (
            (lowland_plain_rank > 0.22)
            | (floodplain_rank > 0.18)
        )
        & (shore_seed | barrier_seed | estuary_seed)
    )
    depositional_domain = (
        low_land
        & (
            delta_plain_mask
            | ((coastal_plain_rank > 0.04) & (distance <= 2.0))
            | supported_lowland
            | process_halo
        )
    )
    depositional_rank = _ranked_halo(
        grid,
        process_seed,
        domain=depositional_domain,
        passes=1,
        decay=0.50,
    )
    depositional_support = np.clip(
        0.34 * coastal_plain_rank
        + 0.24 * lowland_plain_rank
        + 0.22 * floodplain_rank
        + 0.20 * delta_plain_mask.astype(np.float64)
        + 0.14 * np.exp(-np.maximum(distance, 0.0) / 1.55),
        0.0,
        1.0,
    )
    depositional_rank *= np.maximum(depositional_support, 0.36)
    depositional_rank[~depositional_domain] = 0.0
    depositional_rank[depositional_rank < 0.14] = 0.0

    if np.any(barrier_seed & depositional_domain):
        strandplain_rank = _ranked_halo(
            grid,
            barrier_seed & depositional_domain,
            domain=depositional_domain & (distance <= 2.0),
            passes=1,
            decay=0.54,
        )
        strandplain_rank *= np.maximum(depositional_rank, 0.46)
        strandplain_rank[strandplain_rank < 0.12] = 0.0

    tidal_seed = (estuary_seed | delta_plain_mask) & depositional_domain
    if np.any(tidal_seed):
        tidal_flat_rank = _ranked_halo(
            grid,
            tidal_seed,
            domain=depositional_domain & (distance <= 2.2),
            passes=1,
            decay=0.56,
        )
        tidal_flat_rank *= np.maximum(depositional_rank, 0.42)
        tidal_flat_rank[tidal_flat_rank < 0.12] = 0.0

    local_mean = _local_domain_neighbor_mean(grid, refined_rel, land)
    texture = _signed_texture(grid, int(detail_seed) + 6103, scale=6)
    ridged = _ridged_texture(grid, int(detail_seed) + 6119)
    low_gate = np.clip((360.0 - refined_rel) / 360.0, 0.0, 1.0)
    flatten = np.clip(local_mean - refined_rel, -18.0, 18.0)
    plain_delta = depositional_rank * low_gate * (0.34 * flatten + 1.4 * texture)
    strand_delta = strandplain_rank * low_gate * (2.0 + 8.0 * ridged)
    tidal_delta = -tidal_flat_rank * low_gate * (2.5 + 7.5 * (0.55 + 0.45 * ridged))
    delta = np.clip(plain_delta + strand_delta + tidal_delta, -16.0, 18.0)
    delta *= _polar_detail_damping(grid)
    delta[~low_land] = 0.0
    delta[np.abs(delta) < 0.35] = 0.0
    active_deposition = delta != 0.0
    depositional_rank[~active_deposition] = 0.0
    strandplain_rank[~active_deposition] = 0.0
    tidal_flat_rank[~active_deposition] = 0.0
    return (
        delta,
        np.clip(depositional_rank, 0.0, 1.0),
        np.clip(strandplain_rank, 0.0, 1.0),
        np.clip(tidal_flat_rank, 0.0, 1.0),
    )


def _downslope_receivers(
    grid: SphereGrid,
    land: np.ndarray,
    coastal_land: np.ndarray,
    potential: np.ndarray,
    distance: np.ndarray,
) -> np.ndarray:
    receiver = np.full(grid.n, -1, dtype=np.int64)
    receiver[coastal_land] = -2
    for cell in np.flatnonzero(land & ~coastal_land):
        nbs = grid.neighbors[int(cell)]
        land_nbs = nbs[land[nbs]]
        if land_nbs.size == 0:
            receiver[int(cell)] = int(cell)
            continue
        lower_distance = land_nbs[distance[land_nbs] < distance[int(cell)]]
        if lower_distance.size:
            receiver[int(cell)] = int(lower_distance[np.argmin(potential[lower_distance])])
            continue
        best = int(land_nbs[np.argmin(potential[land_nbs])])
        receiver[int(cell)] = best if potential[best] < potential[int(cell)] else int(cell)
    return receiver


def _connect_major_river_paths(
    receiver: np.ndarray,
    land: np.ndarray,
    accumulation: np.ndarray,
    river_rank: np.ndarray,
    *,
    source_threshold: float,
) -> np.ndarray:
    out = np.asarray(river_rank, dtype=np.float64).copy()
    sources = np.flatnonzero(land & (out >= float(source_threshold)))
    if sources.size == 0:
        return out
    sources = sources[np.argsort(accumulation[sources])[::-1]]
    for source in sources[: max(1, min(1200, sources.size))]:
        strength = max(float(out[int(source)]), 0.28)
        cell = int(source)
        seen: set[int] = set()
        for _ in range(512):
            if cell in seen or cell < 0:
                break
            seen.add(cell)
            out[cell] = max(out[cell], strength)
            rec = int(receiver[cell])
            if rec < 0 or rec == cell:
                break
            cell = rec
            strength = max(0.24, strength * 0.993)
    return out


def _major_river_path_rank(
    grid: SphereGrid,
    receiver: np.ndarray,
    land: np.ndarray,
    coastal_land: np.ndarray,
    accumulation: np.ndarray,
    distance: np.ndarray,
    river_rank: np.ndarray,
) -> np.ndarray:
    """Extract sparse basin-scale river objects from the dense flow field."""
    out = np.zeros(grid.n, dtype=np.float64)
    land_cells = np.flatnonzero(land)
    if land_cells.size == 0:
        return out

    finite_distance = distance[land_cells]
    d95 = max(float(np.nanpercentile(finite_distance, 95.0)), 1.0)
    min_source_distance = max(3.0, min(7.0, 0.28 * d95))
    acc_land = np.asarray(accumulation[land_cells], dtype=np.float64)
    acc_lo = float(np.nanpercentile(acc_land, 82.0))
    acc_hi = max(float(np.nanpercentile(acc_land, 99.5)), acc_lo + 1.0e-9)
    acc_rank = np.zeros(grid.n, dtype=np.float64)
    acc_rank[land_cells] = np.clip(
        (np.asarray(accumulation[land_cells], dtype=np.float64) - acc_lo)
        / (acc_hi - acc_lo),
        0.0,
        1.0,
    )
    inland_bonus = np.clip(distance / max(d95, 1.0), 0.0, 1.0) ** 0.45
    source_score = (0.55 * np.clip(river_rank, 0.0, 1.0) + 0.45 * acc_rank)
    source_score *= 0.35 + 0.65 * inland_bonus
    candidates = (
        land
        & ~coastal_land
        & (distance >= min_source_distance)
        & (source_score >= 0.18)
    )
    sources = np.flatnonzero(candidates)
    if sources.size == 0:
        return out
    sources = sources[np.argsort(source_score[sources])[::-1]]

    max_paths = int(np.clip(round(0.0030 * land_cells.size), 12, 110))
    min_path_len = int(np.clip(round(0.24 * d95), 3, 10))
    blocked_sources = np.zeros(grid.n, dtype=bool)
    blocked_mouths = np.zeros(grid.n, dtype=bool)
    selected = 0
    for source in sources:
        source = int(source)
        if selected >= max_paths:
            break
        if blocked_sources[source]:
            continue
        path = _trace_receiver_path(receiver, source, max_steps=512)
        if len(path) < min_path_len:
            continue
        mouth = int(path[-1])
        if blocked_mouths[mouth]:
            continue

        strength = max(float(np.nanmax(river_rank[path])), float(source_score[source]), 0.42)
        denom = max(len(path) - 1, 1)
        for i, cell in enumerate(path):
            downstream = i / denom
            out[int(cell)] = max(out[int(cell)], strength * (0.72 + 0.28 * downstream))

        path_mask = np.zeros(grid.n, dtype=bool)
        path_mask[np.asarray(path, dtype=np.int64)] = True
        blocked_sources |= _dilate_mask(grid, path_mask, passes=2)
        mouth_mask = np.zeros(grid.n, dtype=bool)
        mouth_mask[mouth] = True
        blocked_mouths |= _dilate_mask(grid, mouth_mask, passes=3)
        selected += 1
    if selected == 0:
        return np.clip(out, 0.0, 1.0)

    # Long basin-outlet A* rivers are intentionally not part of the default
    # selected-snapshot pass yet.  Early trials made longer rivers, but also
    # introduced looped edge paths and a 3x runtime hit.  Keep the default
    # conservative until drainage basins are explicit objects.
    main_mask = out > 0.0
    tributary_candidates = np.flatnonzero(
        land
        & ~coastal_land
        & ~main_mask
        & (distance >= max(2.0, 0.18 * d95))
        & (river_rank >= 0.40)
    )
    if tributary_candidates.size:
        tributary_candidates = tributary_candidates[
            np.argsort(source_score[tributary_candidates])[::-1]
        ]
    selected_tributaries = 0
    max_tributaries = max_paths * 3
    for source in tributary_candidates:
        source = int(source)
        if selected_tributaries >= max_tributaries:
            break
        if blocked_sources[source]:
            continue
        path = _trace_receiver_path(receiver, source, max_steps=256)
        if len(path) < 3:
            continue
        join_index = -1
        for i, cell in enumerate(path[1:], start=1):
            if main_mask[int(cell)]:
                join_index = i
                break
        if join_index < 2:
            continue
        tributary = path[:join_index]
        strength = min(0.66, max(0.30, float(np.nanmax(river_rank[tributary])) * 0.72))
        denom = max(len(tributary) - 1, 1)
        for i, cell in enumerate(tributary):
            downstream = i / denom
            out[int(cell)] = max(out[int(cell)], strength * (0.70 + 0.30 * downstream))
        path_mask = np.zeros(grid.n, dtype=bool)
        path_mask[np.asarray(tributary, dtype=np.int64)] = True
        blocked_sources |= _dilate_mask(grid, path_mask, passes=1)
        selected_tributaries += 1
    return np.clip(out, 0.0, 1.0)


def _drainage_basin_objects(
    grid: SphereGrid,
    receiver: np.ndarray,
    land: np.ndarray,
    coastal_land: np.ndarray,
    accumulation: np.ndarray,
    distance: np.ndarray,
    river_rank: np.ndarray,
    rel: np.ndarray,
    basin_region: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Select explicit basin trunk/floodplain objects from the receiver tree."""
    basin_id = np.zeros(grid.n, dtype=np.int32)
    trunk_rank = np.zeros(grid.n, dtype=np.float64)
    floodplain_rank = np.zeros(grid.n, dtype=np.float64)
    land_cells = np.flatnonzero(land)
    if land_cells.size == 0:
        return basin_id, trunk_rank, floodplain_rank

    finite_distance = np.asarray(distance[land_cells], dtype=np.float64)
    d95 = max(float(np.nanpercentile(finite_distance, 95.0)), 1.0)
    dmax = max(float(np.nanmax(finite_distance)), d95)
    acc_land = np.asarray(accumulation[land_cells], dtype=np.float64)
    acc_lo = float(np.nanpercentile(acc_land, 78.0))
    acc_hi = max(float(np.nanpercentile(acc_land, 99.5)), acc_lo + 1.0e-9)
    acc_rank = np.zeros(grid.n, dtype=np.float64)
    acc_rank[land_cells] = np.clip(
        (np.asarray(accumulation[land_cells], dtype=np.float64) - acc_lo)
        / (acc_hi - acc_lo),
        0.0,
        1.0,
    )
    dist_rank = np.zeros(grid.n, dtype=np.float64)
    dist_rank[land_cells] = np.clip(distance[land_cells] / max(d95, 1.0),
                                    0.0, 1.0)
    lowland_rank = np.zeros(grid.n, dtype=np.float64)
    lowland_rank[land_cells] = np.clip((950.0 - rel[land_cells]) / 950.0,
                                       0.0,
                                       1.0)
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    basin_lat_weight = np.clip((70.0 - lat_abs) / 12.0, 0.0, 1.0)
    candidate_score = (
        0.36 * acc_rank
        + 0.34 * dist_rank
        + 0.18 * np.clip(river_rank, 0.0, 1.0)
        + 0.08 * lowland_rank
        + 0.04 * np.asarray(basin_region, dtype=np.float64)
    )
    candidate_score *= basin_lat_weight
    candidates = np.flatnonzero(
        land
        & ~coastal_land
        & (lat_abs < 66.0)
        & (distance >= max(4.0, 0.36 * d95))
        & (candidate_score >= 0.18)
    )
    if candidates.size == 0:
        return basin_id, trunk_rank, floodplain_rank

    candidates = candidates[np.argsort(candidate_score[candidates])[::-1]]
    max_trace = min(int(candidates.size), 5200)
    min_path_len = int(np.clip(round(0.48 * d95), 5, 10))
    path_records: list[tuple[float, int, int, list[int], float]] = []
    for source in candidates[:max_trace]:
        source = int(source)
        path = _trace_receiver_path(receiver, source, max_steps=512)
        if len(path) < min_path_len:
            continue
        mouth = int(path[-1])
        if mouth < 0 or not coastal_land[mouth]:
            continue
        path_arr = np.asarray(path, dtype=np.int64)
        if abs(float(grid.lat[mouth])) >= 70.0:
            continue
        if float(np.nanpercentile(lat_abs[path_arr], 90.0)) >= 68.0:
            continue
        path_len_rank = np.clip((len(path) - min_path_len)
                                / max(dmax - min_path_len, 1.0),
                                0.0,
                                1.0)
        path_dist_rank = np.clip(np.nanmax(distance[path_arr]) / max(d95, 1.0),
                                 0.0,
                                 1.0)
        mouth_rank = np.clip(acc_rank[mouth], 0.0, 1.0)
        lowland_path = float(np.nanmean(lowland_rank[path_arr]))
        priority = (
            0.36 * float(candidate_score[source])
            + 0.30 * float(path_len_rank)
            + 0.18 * float(path_dist_rank)
            + 0.10 * float(mouth_rank)
            + 0.06 * lowland_path
        )
        strength = max(0.42, min(1.0, 0.50 + 0.58 * priority))
        path_records.append((priority, source, mouth, path, strength))
    if not path_records:
        return basin_id, trunk_rank, floodplain_rank
    path_records.sort(key=lambda item: item[0], reverse=True)

    max_basins = int(np.clip(round(0.00085 * land_cells.size), 10, 28))
    blocked_sources = np.zeros(grid.n, dtype=bool)
    blocked_mouths = np.zeros(grid.n, dtype=bool)
    selected_paths: list[tuple[int, int, list[int], float]] = []
    selected_mouths: dict[int, int] = {}
    for _, source, mouth, path, strength in path_records:
        if len(selected_paths) >= max_basins:
            break
        if blocked_sources[source] or blocked_mouths[mouth]:
            continue
        path_arr = np.asarray(path, dtype=np.int64)
        if np.count_nonzero(trunk_rank[path_arr] > 0.0) > max(2, int(0.30 * len(path))):
            continue
        basin_num = len(selected_paths) + 1
        selected_paths.append((basin_num, mouth, path, strength))
        selected_mouths[mouth] = basin_num
        denom = max(len(path) - 1, 1)
        for i, cell in enumerate(path):
            downstream = i / denom
            trunk_rank[int(cell)] = max(trunk_rank[int(cell)],
                                        strength * (0.68 + 0.32 * downstream))
        path_mask = np.zeros(grid.n, dtype=bool)
        path_mask[path_arr] = True
        blocked_sources |= _dilate_mask(grid, path_mask, passes=3)
        mouth_mask = np.zeros(grid.n, dtype=bool)
        mouth_mask[mouth] = True
        blocked_mouths |= _dilate_mask(grid, mouth_mask, passes=5)

    if not selected_paths:
        return basin_id, trunk_rank, floodplain_rank

    mouth_cache: dict[int, int] = {}
    for cell in land_cells:
        cell = int(cell)
        mouth = _receiver_mouth(receiver, cell, mouth_cache)
        if mouth in selected_mouths:
            basin_id[cell] = selected_mouths[mouth]

    lowland = land & (lat_abs < 68.0) & ((rel < 900.0) | basin_region) & (distance >= 1.0)
    for basin_num, _, path, strength in selected_paths:
        path_arr = np.asarray(path, dtype=np.int64)
        path_mask = np.zeros(grid.n, dtype=bool)
        path_mask[path_arr] = True
        near = _dilate_mask(grid, path_mask, passes=1) & lowland
        broad = _dilate_mask(grid, path_mask, passes=2) & lowland & (rel < 650.0)
        floodplain_rank[near] = np.maximum(floodplain_rank[near], 0.46 * strength)
        floodplain_rank[broad] = np.maximum(floodplain_rank[broad], 0.24 * strength)
        basin_id[(near | broad) & (basin_id == 0)] = basin_num
    floodplain_rank[trunk_rank > 0.0] = np.maximum(
        floodplain_rank[trunk_rank > 0.0],
        0.62 * trunk_rank[trunk_rank > 0.0],
    )
    return basin_id, np.clip(trunk_rank, 0.0, 1.0), np.clip(floodplain_rank, 0.0, 1.0)


def _meander_belt_rank(
    grid: SphereGrid,
    land: np.ndarray,
    rel: np.ndarray,
    distance: np.ndarray,
    basin_trunk_rank: np.ndarray,
    floodplain_rank: np.ndarray,
    basin_region: np.ndarray,
) -> np.ndarray:
    """Proxy for low-gradient meander belts around selected trunk rivers."""
    trunk = np.asarray(basin_trunk_rank, dtype=np.float64) > 0.0
    if not trunk.any():
        return np.zeros(grid.n, dtype=np.float64)
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    low_polar = np.clip((66.0 - lat_abs) / 14.0, 0.0, 1.0)
    lowland = np.clip((720.0 - np.asarray(rel, dtype=np.float64)) / 720.0,
                      0.0,
                      1.0)
    inland = np.clip(np.asarray(distance, dtype=np.float64) / 8.0, 0.0, 1.0)
    base = (
        0.58 * np.asarray(floodplain_rank, dtype=np.float64)
        + 0.42 * np.asarray(basin_trunk_rank, dtype=np.float64)
    )
    belt = base * lowland * (0.45 + 0.55 * inland) * low_polar
    belt *= np.where(np.asarray(basin_region, dtype=bool), 1.12, 1.0)
    belt[~land] = 0.0
    belt[rel >= 950.0] = 0.0
    belt[lat_abs >= 68.0] = 0.0
    return np.clip(belt, 0.0, 1.0)


def _receiver_mouth(
    receiver: np.ndarray,
    source: int,
    cache: dict[int, int],
) -> int:
    source = int(source)
    if source in cache:
        return cache[source]
    path: list[int] = []
    cell = source
    seen: set[int] = set()
    mouth = -1
    for _ in range(512):
        if cell in cache:
            mouth = int(cache[cell])
            break
        if cell < 0 or cell in seen:
            mouth = -1
            break
        path.append(cell)
        seen.add(cell)
        rec = int(receiver[cell])
        if rec < 0:
            mouth = cell
            break
        if rec == cell:
            mouth = -1
            break
        cell = rec
    for item in path:
        cache[int(item)] = int(mouth)
    return int(mouth)


def _trace_receiver_path(
    receiver: np.ndarray,
    source: int,
    *,
    max_steps: int,
) -> list[int]:
    path: list[int] = []
    cell = int(source)
    seen: set[int] = set()
    for _ in range(max(1, int(max_steps))):
        if cell < 0 or cell in seen:
            break
        path.append(cell)
        seen.add(cell)
        rec = int(receiver[cell])
        if rec < 0 or rec == cell:
            break
        cell = rec
    return path


def _land_distance_from_ocean(grid: SphereGrid, land: np.ndarray) -> np.ndarray:
    distance = np.full(grid.n, np.inf, dtype=np.float64)
    queue: deque[int] = deque()
    for cell in np.flatnonzero(land):
        nbs = grid.neighbors[int(cell)]
        if nbs.size and np.any(~land[nbs]):
            distance[int(cell)] = 0.0
            queue.append(int(cell))
    while queue:
        cell = queue.popleft()
        next_distance = distance[cell] + 1.0
        for nb in grid.neighbors[cell]:
            nb = int(nb)
            if land[nb] and next_distance < distance[nb]:
                distance[nb] = next_distance
                queue.append(nb)
    return distance


def _ocean_distance_from_land(grid: SphereGrid, ocean: np.ndarray) -> np.ndarray:
    distance = np.full(grid.n, np.inf, dtype=np.float64)
    queue: deque[int] = deque()
    land = ~np.asarray(ocean, dtype=bool)
    for cell in np.flatnonzero(ocean):
        nbs = grid.neighbors[int(cell)]
        if nbs.size and np.any(land[nbs]):
            distance[int(cell)] = 0.0
            queue.append(int(cell))
    while queue:
        cell = queue.popleft()
        next_distance = distance[cell] + 1.0
        for nb in grid.neighbors[cell]:
            nb = int(nb)
            if ocean[nb] and next_distance < distance[nb]:
                distance[nb] = next_distance
                queue.append(nb)
    return distance


def _ocean_coast_distance(
    grid: SphereGrid,
    land: np.ndarray,
    ocean: np.ndarray,
    *,
    max_passes: int,
) -> np.ndarray:
    distance = np.full(grid.n, -1, dtype=np.int16)
    current = _ocean_neighbors_of(grid, land, ocean)
    distance[current] = 0
    for step in range(1, int(max_passes) + 1):
        grown = _dilate_mask(grid, current, passes=1) & ocean & (distance < 0)
        if not grown.any():
            break
        distance[grown] = step
        current |= grown
    return distance


def _ocean_neighbors_of(
    grid: SphereGrid,
    source: np.ndarray,
    ocean: np.ndarray,
) -> np.ndarray:
    out = np.zeros(grid.n, dtype=bool)
    for cell in np.flatnonzero(source):
        nbs = grid.neighbors[int(cell)]
        out[nbs[ocean[nbs]]] = True
    return out


def _land_neighbors_of(
    grid: SphereGrid,
    source: np.ndarray,
    land: np.ndarray,
) -> np.ndarray:
    out = np.zeros(grid.n, dtype=bool)
    for cell in np.flatnonzero(source):
        nbs = grid.neighbors[int(cell)]
        out[nbs[land[nbs]]] = True
    return out


def _open_ocean_candidate_weight(
    ocean_distance: np.ndarray,
    *,
    start: float,
    full: float,
) -> np.ndarray:
    d = np.asarray(ocean_distance, dtype=np.float64)
    span = max(float(full) - float(start), 1.0e-9)
    return np.clip((d - float(start)) / span, 0.0, 1.0)


def _offshore_island_candidate_gate(
    ocean_distance: np.ndarray,
    *,
    seamount: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    island_arc: np.ndarray,
) -> np.ndarray:
    d = np.asarray(ocean_distance, dtype=np.float64)
    return (
        (np.asarray(seamount, dtype=bool) & (d >= 2.0))
        | (np.asarray(island_arc, dtype=bool) & (d >= 2.0))
        | (np.asarray(plateau, dtype=bool) & (d >= 3.0))
        | (np.asarray(microcontinent, dtype=bool) & (d >= 4.0))
    )


def _offshore_atoll_candidate_gate(
    ocean_distance: np.ndarray,
    *,
    seamount: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    island_arc: np.ndarray,
) -> np.ndarray:
    d = np.asarray(ocean_distance, dtype=np.float64)
    return (
        (np.asarray(seamount, dtype=bool) & (d >= 3.0))
        | (np.asarray(island_arc, dtype=bool) & (d >= 4.0))
        | (np.asarray(plateau, dtype=bool) & (d >= 4.0))
        | (np.asarray(microcontinent, dtype=bool) & (d >= 5.0))
    )


def _marine_island_candidates(
    grid: SphereGrid,
    rel: np.ndarray,
    potential_rel: np.ndarray,
    process_high: np.ndarray,
    trench: np.ndarray,
    seamount: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    island_arc: np.ndarray,
    ridge_texture: np.ndarray,
    low_polar: np.ndarray,
    ocean_distance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Rank process-backed islet candidates without changing ocean topology."""
    rel = np.asarray(rel, dtype=np.float64)
    potential_rel = np.asarray(potential_rel, dtype=np.float64)
    ocean = rel < 0.0
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    lat_suitability = np.clip((62.0 - lat_abs) / 22.0, 0.0, 1.0)
    d = np.asarray(ocean_distance, dtype=np.float64)
    offshore_gate = _offshore_island_candidate_gate(
        d,
        seamount=seamount,
        plateau=plateau,
        microcontinent=microcontinent,
        island_arc=island_arc,
    )
    domain = (
        ocean
        & np.asarray(process_high, dtype=bool)
        & np.asarray(low_polar, dtype=bool)
        & ~np.asarray(trench, dtype=bool)
        & offshore_gate
        & (rel > -1750.0)
        & (potential_rel > -820.0)
        & (lat_abs <= 62.0)
    )
    rank = np.zeros(grid.n, dtype=np.float64)
    mask = np.zeros(grid.n, dtype=bool)
    if not domain.any():
        return rank, mask

    depth = np.maximum(-rel, 0.0)
    shallowness = np.clip((1500.0 - depth) / 1500.0, 0.0, 1.0)
    emergence = np.clip((potential_rel + 520.0) / 520.0, 0.0, 1.0)
    uplift = np.clip((potential_rel - rel) / 460.0, 0.0, 1.0)
    open_water = _open_ocean_candidate_weight(d, start=2.0, full=8.0)
    process_bonus = (
        0.06 * np.asarray(seamount, dtype=np.float64)
        + 0.10 * np.asarray(plateau, dtype=np.float64)
        + 0.15 * np.asarray(microcontinent, dtype=np.float64)
        + 0.13 * np.asarray(island_arc, dtype=np.float64)
    )
    score = (
        0.38 * emergence
        + 0.26 * shallowness
        + 0.20 * uplift
        + 0.12 * np.asarray(ridge_texture, dtype=np.float64)
        + 0.12 * open_water
        + process_bonus
    )
    score *= 0.55 + 0.45 * lat_suitability
    score *= 0.62 + 0.38 * open_water
    score = np.clip(score, 0.0, 1.0)
    eligible = domain & (score >= 0.50)
    if not eligible.any():
        eligible = domain & (score >= max(0.42, float(np.nanpercentile(score[domain], 88.0))))
    if not eligible.any():
        return rank, mask

    count = min(
        max(6, int(round(0.00045 * grid.n))),
        int(np.count_nonzero(eligible)),
    )
    chosen = np.flatnonzero(eligible)[np.argsort(score[eligible])[-count:]]
    mask[chosen] = True
    halo_domain = domain & (score >= 0.36)
    rank = _ranked_halo(grid, mask, domain=halo_domain, passes=1, decay=0.54)
    rank *= np.maximum(score, 0.62)
    rank[~halo_domain] = 0.0
    rank[mask] = np.maximum(rank[mask], score[mask])
    return np.clip(rank, 0.0, 1.0), mask


def _marine_reef_atoll_morphology(
    grid: SphereGrid,
    rel: np.ndarray,
    potential_rel: np.ndarray,
    reef_domain: np.ndarray,
    reef_atoll: np.ndarray,
    atoll_candidate: np.ndarray,
    atoll_candidate_rank: np.ndarray,
    process_island_candidate: np.ndarray,
    process_high: np.ndarray,
    microcontinent: np.ndarray,
    island_arc: np.ndarray,
    plateau: np.ndarray,
    trench: np.ndarray,
    ocean_distance: np.ndarray,
    ridge_texture: np.ndarray,
    tropical: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Derive local reef rim, lagoon, and fringing-reef ranks."""
    rel = np.asarray(rel, dtype=np.float64)
    potential_rel = np.asarray(potential_rel, dtype=np.float64)
    ocean = rel < 0.0
    reef_domain = np.asarray(reef_domain, dtype=bool)
    atoll_candidate = np.asarray(atoll_candidate, dtype=bool)
    atoll_candidate_rank = np.asarray(atoll_candidate_rank, dtype=np.float64)
    process_high = np.asarray(process_high, dtype=bool)
    trench = np.asarray(trench, dtype=bool)
    d = np.asarray(ocean_distance, dtype=np.float64)
    ridge_texture = np.asarray(ridge_texture, dtype=np.float64)

    lagoon_rank = np.zeros(grid.n, dtype=np.float64)
    reef_rim_rank = np.zeros(grid.n, dtype=np.float64)
    fringing_rank = np.zeros(grid.n, dtype=np.float64)

    atoll_values = atoll_candidate_rank[atoll_candidate]
    if atoll_values.size:
        core_cut = max(0.72, float(np.nanpercentile(atoll_values, 82.0)))
    else:
        core_cut = 1.1
    lagoon_seed = atoll_candidate & (atoll_candidate_rank >= core_cut)
    if lagoon_seed.any():
        lagoon_domain = (
            reef_domain
            & ocean
            & ~trench
            & (rel > -950.0)
            & (potential_rel > -420.0)
        )
        lagoon_rank = _ranked_halo(
            grid,
            lagoon_seed,
            domain=lagoon_domain,
            passes=1,
            decay=0.34,
        )
        lagoon_rank *= np.maximum(atoll_candidate_rank, 0.52)
        lagoon_rank[~lagoon_domain] = 0.0
        rim_seed = (
            (atoll_candidate & ~lagoon_seed)
            | (_dilate_mask(grid, lagoon_seed, passes=1) & lagoon_domain & ~lagoon_seed)
        )
        reef_rim_rank = _ranked_halo(
            grid,
            rim_seed,
            domain=lagoon_domain,
            passes=1,
            decay=0.55,
        )
        reef_rim_rank *= np.maximum(0.58, atoll_candidate_rank)
        reef_rim_rank[lagoon_rank >= 0.55] *= 0.18

    shallow = np.clip((620.0 + rel) / 620.0, 0.0, 1.0)
    nearshore = np.clip((2.35 - d) / 2.35, 0.0, 1.0)
    process_bonus = (
        0.16 * np.asarray(microcontinent, dtype=np.float64)
        + 0.13 * np.asarray(island_arc, dtype=np.float64)
        + 0.08 * np.asarray(plateau, dtype=np.float64)
        + 0.10 * np.asarray(process_island_candidate, dtype=np.float64)
        + 0.08 * np.asarray(reef_atoll, dtype=np.float64)
    )
    fringing_score = (
        0.38 * shallow
        + 0.28 * nearshore
        + 0.16 * ridge_texture
        + process_bonus
    )
    fringing_domain = (
        ocean
        & np.asarray(tropical, dtype=bool)
        & ~trench
        & process_high
        & (d >= 0.0)
        & (d <= 2.35)
        & (rel > -680.0)
        & (rel < -4.0)
        & (fringing_score >= 0.52)
    )
    if fringing_domain.any():
        fringing_seed = fringing_domain & (fringing_score >= max(
            0.58,
            float(np.nanpercentile(fringing_score[fringing_domain], 70.0)),
        ))
        fringing_rank = _ranked_halo(
            grid,
            fringing_seed,
            domain=fringing_domain,
            passes=1,
            decay=0.50,
        )
        fringing_rank *= np.maximum(fringing_score, 0.58)
        fringing_rank[~fringing_domain] = 0.0

    return (
        np.clip(reef_rim_rank, 0.0, 1.0),
        np.clip(lagoon_rank, 0.0, 1.0),
        np.clip(fringing_rank, 0.0, 1.0),
    )


def _reef_atoll_morphology_delta(
    rel: np.ndarray,
    reef_rim_rank: np.ndarray,
    atoll_lagoon_rank: np.ndarray,
    fringing_reef_rank: np.ndarray,
) -> np.ndarray:
    """Apply small bathymetric relief for reefs while preserving ocean sign."""
    rel = np.asarray(rel, dtype=np.float64)
    delta = np.zeros(rel.shape, dtype=np.float64)

    rim = np.asarray(reef_rim_rank, dtype=np.float64)
    rim_target = -18.0 - 8.0 * (1.0 - rim)
    rim_lift = rim > 0.05
    delta[rim_lift] += np.minimum(
        54.0,
        np.maximum(0.0, rim_target[rim_lift] - rel[rim_lift]) * 0.52,
    )

    fringing = np.asarray(fringing_reef_rank, dtype=np.float64)
    fringing_target = -22.0 - 10.0 * (1.0 - fringing)
    fringing_lift = fringing > 0.05
    delta[fringing_lift] += np.minimum(
        42.0,
        np.maximum(0.0, fringing_target[fringing_lift] - rel[fringing_lift]) * 0.42,
    )

    lagoon = np.asarray(atoll_lagoon_rank, dtype=np.float64)
    lagoon_target = -56.0 - 34.0 * (1.0 - lagoon)
    lagoon_deepen = lagoon > 0.05
    delta[lagoon_deepen] -= np.minimum(
        36.0,
        np.maximum(0.0, rel[lagoon_deepen] - lagoon_target[lagoon_deepen]) * 0.38,
    )
    return delta


def _marine_process_island_promotion(
    grid: SphereGrid,
    rel: np.ndarray,
    potential_rel: np.ndarray,
    island_candidate_rank: np.ndarray,
    atoll_candidate_rank: np.ndarray,
    reef_rim_rank: np.ndarray,
    atoll_lagoon_rank: np.ndarray,
    fringing_reef_rank: np.ndarray,
    process_island_candidate: np.ndarray,
    atoll_candidate: np.ndarray,
    seamount: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    island_arc: np.ndarray,
    trench: np.ndarray,
    low_polar: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Promote a very small number of process-backed islets when explicitly enabled."""
    rel = np.asarray(rel, dtype=np.float64)
    potential_rel = np.asarray(potential_rel, dtype=np.float64)
    ocean = rel < 0.0
    island_rank = np.asarray(island_candidate_rank, dtype=np.float64)
    atoll_rank = np.asarray(atoll_candidate_rank, dtype=np.float64)
    rim_rank = np.asarray(reef_rim_rank, dtype=np.float64)
    lagoon_rank = np.asarray(atoll_lagoon_rank, dtype=np.float64)
    fringing_rank = np.asarray(fringing_reef_rank, dtype=np.float64)
    process_high = (
        np.asarray(seamount, dtype=bool)
        | np.asarray(plateau, dtype=bool)
        | np.asarray(microcontinent, dtype=bool)
        | np.asarray(island_arc, dtype=bool)
    )
    trench = np.asarray(trench, dtype=bool)
    low_polar = np.asarray(low_polar, dtype=bool)

    delta = np.zeros(grid.n, dtype=np.float64)
    promotion_rank = np.zeros(grid.n, dtype=np.float64)
    promoted = np.zeros(grid.n, dtype=bool)
    atoll_promoted = np.zeros(grid.n, dtype=bool)

    reef_land_rank = np.maximum(rim_rank, 0.78 * fringing_rank)
    atoll_land_score = reef_land_rank * np.maximum(atoll_rank, 0.52)
    score = np.maximum(island_rank, atoll_land_score)
    base_domain = (
        ocean
        & low_polar
        & ~trench
        & process_high
        & (potential_rel > -160.0)
    )
    islet_domain = (
        base_domain
        & np.asarray(process_island_candidate, dtype=bool)
        & (island_rank >= 0.66)
    )
    atoll_domain = (
        base_domain
        & np.asarray(atoll_candidate, dtype=bool)
        & (atoll_rank >= 0.60)
        & (reef_land_rank >= 0.34)
        & (lagoon_rank < 0.80)
        & (potential_rel > -120.0)
    )
    candidate_domain = islet_domain | atoll_domain
    if not candidate_domain.any():
        return delta, promotion_rank, promoted, atoll_promoted

    seed = np.zeros(grid.n, dtype=bool)
    if islet_domain.any():
        islet_cut = max(0.69, float(np.nanpercentile(island_rank[islet_domain], 78.0)))
        islet_eligible = islet_domain & (island_rank >= islet_cut)
        islet_count = min(
            max(2, int(round(0.00013 * grid.n))),
            int(np.count_nonzero(islet_eligible)),
        )
        if islet_count > 0:
            cells = np.flatnonzero(islet_eligible)
            seed[cells[np.argsort(island_rank[islet_eligible])[-islet_count:]]] = True
    if atoll_domain.any():
        atoll_score = np.maximum(atoll_land_score, 0.58 * atoll_rank)
        atoll_cut = max(0.50, float(np.nanpercentile(atoll_score[atoll_domain], 68.0)))
        atoll_eligible = atoll_domain & (atoll_score >= atoll_cut)
        atoll_count = min(
            max(1, int(round(0.000045 * grid.n))),
            int(np.count_nonzero(atoll_eligible)),
        )
        if atoll_count > 0:
            cells = np.flatnonzero(atoll_eligible)
            seed[cells[np.argsort(atoll_score[atoll_eligible])[-atoll_count:]]] = True
    if not seed.any():
        return delta, promotion_rank, promoted, atoll_promoted

    halo = (
        _dilate_mask(grid, seed, passes=1)
        & candidate_domain
        & (
            (islet_domain & (island_rank >= 0.73))
            | (atoll_domain & (atoll_land_score >= 0.44))
        )
    )
    promoted = seed | halo
    max_promoted = min(
        max(4, int(round(0.00028 * grid.n))),
        int(np.count_nonzero(candidate_domain)),
    )
    if np.count_nonzero(promoted) > max_promoted:
        cells = np.flatnonzero(promoted)
        keep = cells[np.argsort(score[promoted])[-max_promoted:]]
        limited = np.zeros(grid.n, dtype=bool)
        limited[keep] = True
        promoted = limited

    promotion_rank = _ranked_halo(
        grid,
        promoted,
        domain=candidate_domain,
        passes=1,
        decay=0.45,
    )
    promotion_rank *= np.maximum(score, 0.58)
    promotion_rank[promoted] = np.maximum(promotion_rank[promoted], score[promoted])
    promotion_rank[~candidate_domain] = 0.0
    atoll_promoted = (
        promoted
        & np.asarray(atoll_candidate, dtype=bool)
        & (reef_land_rank >= np.maximum(0.34, 0.72 * island_rank))
        & (lagoon_rank < 0.82)
    )

    oceanic_island_target = (
        18.0
        + 105.0 * island_rank
        + 38.0 * np.asarray(microcontinent, dtype=np.float64)
        + 26.0 * np.asarray(island_arc, dtype=np.float64)
        + 14.0 * np.asarray(plateau, dtype=np.float64)
    )
    reef_islet_target = 4.0 + 28.0 * np.maximum(rim_rank, fringing_rank)
    target = np.where(atoll_promoted, reef_islet_target, oceanic_island_target)
    target = np.clip(target, 3.0, 190.0)
    delta[promoted] = np.maximum(0.0, target[promoted] - potential_rel[promoted])
    return (
        delta,
        np.clip(promotion_rank, 0.0, 1.0),
        promoted,
        atoll_promoted,
    )


def _marine_island_microshape_ranks(
    island_candidate_rank: np.ndarray,
    atoll_candidate_rank: np.ndarray,
    reef_rim_rank: np.ndarray,
    atoll_lagoon_rank: np.ndarray,
    fringing_reef_rank: np.ndarray,
    process_island_candidate: np.ndarray,
    atoll_candidate: np.ndarray,
    process_island_promotion_rank: np.ndarray,
    process_island_promoted: np.ndarray,
    atoll_islet_promoted: np.ndarray,
    ocean: np.ndarray,
    low_polar: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Ranks used only by subcell visual QA for islands and atolls."""
    ocean = np.asarray(ocean, dtype=bool)
    low_polar = np.asarray(low_polar, dtype=bool)
    island_rank = np.asarray(island_candidate_rank, dtype=np.float64)
    atoll_rank = np.asarray(atoll_candidate_rank, dtype=np.float64)
    promotion_rank = np.asarray(process_island_promotion_rank, dtype=np.float64)
    reef_land = np.maximum(
        np.asarray(reef_rim_rank, dtype=np.float64),
        0.78 * np.asarray(fringing_reef_rank, dtype=np.float64),
    )
    lagoon = np.asarray(atoll_lagoon_rank, dtype=np.float64)

    islet = np.zeros(island_rank.shape, dtype=np.float64)
    islet_domain = (
        ocean
        & low_polar
        & (
            np.asarray(process_island_candidate, dtype=bool)
            | np.asarray(process_island_promoted, dtype=bool)
        )
    )
    islet[islet_domain] = np.maximum(
        island_rank[islet_domain],
        promotion_rank[islet_domain],
    )

    atoll = np.zeros(atoll_rank.shape, dtype=np.float64)
    atoll_domain = (
        ocean
        & low_polar
        & (
            np.asarray(atoll_candidate, dtype=bool)
            | np.asarray(atoll_islet_promoted, dtype=bool)
        )
    )
    atoll_score = np.maximum(
        atoll_rank * np.maximum(reef_land, 0.42),
        promotion_rank * np.asarray(atoll_islet_promoted, dtype=np.float64),
    )
    atoll_score *= np.where(lagoon >= 0.80, 0.64, 1.0)
    atoll[atoll_domain] = atoll_score[atoll_domain]
    return np.clip(islet, 0.0, 1.0), np.clip(atoll, 0.0, 1.0)


def _marine_shelf_slope_microrelief_delta(
    grid: SphereGrid,
    rel: np.ndarray,
    refined_rel: np.ndarray,
    shelf_break_rank: np.ndarray,
    ocean_distance: np.ndarray,
    passive_margin: np.ndarray,
    active_margin: np.ndarray,
    trench: np.ndarray,
    ridge: np.ndarray,
    process_island_promoted: np.ndarray,
    low_polar: np.ndarray,
    detail_seed: int,
) -> np.ndarray:
    """Add sparse shelf-edge and upper-slope microrelief without changing topology."""
    rel = np.asarray(rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    ocean = rel < 0.0
    shelf_break = np.asarray(shelf_break_rank, dtype=np.float64)
    d = np.asarray(ocean_distance, dtype=np.float64)
    passive_margin = np.asarray(passive_margin, dtype=bool)
    active_margin = np.asarray(active_margin, dtype=bool)
    trench = np.asarray(trench, dtype=bool)
    ridge = np.asarray(ridge, dtype=bool)
    low_polar = np.asarray(low_polar, dtype=bool)
    promoted = np.asarray(process_island_promoted, dtype=bool)
    domain = (
        ocean
        & low_polar
        & ~trench
        & ~ridge
        & ~promoted
        & (d >= 0.0)
        & (d <= 7.0)
        & (shelf_break > 0.42)
        & (rel > -2800.0)
        & (refined_rel < -2.0)
    )
    delta = np.zeros(rel.shape, dtype=np.float64)
    if not np.any(domain):
        return delta

    texture = _signed_texture(grid, int(detail_seed) + 3533, scale=5)
    ridged = _ridged_texture(grid, int(detail_seed) + 3547)
    depth_gate = np.clip((3100.0 + refined_rel) / 3100.0, 0.0, 1.0)
    margin_support = (
        0.72
        + 0.18 * passive_margin.astype(np.float64)
        + 0.12 * active_margin.astype(np.float64)
    )
    score = np.clip(
        shelf_break
        * margin_support
        * (0.62 + 0.38 * ridged)
        * depth_gate,
        0.0,
        1.0,
    )
    score[~domain] = 0.0
    rank = _shelf_slope_axis_rank(grid, domain, score, shelf_break)
    if not np.any(rank > 0.0):
        return delta

    channel_seed = (
        domain
        & (rank >= 0.30)
        & (texture < -0.10)
        & (score >= 0.48)
    )
    channel_support = domain & (rank >= 0.12) & (texture < 0.12)
    channel_axis = _bridge_sparse_line_mask(
        grid,
        channel_seed,
        domain=channel_support,
        passes=1,
    )
    channel_rank = np.zeros(rel.shape, dtype=np.float64)
    if np.any(channel_axis):
        channel_rank = _ranked_halo(
            grid,
            channel_axis,
            domain=channel_support,
            passes=1,
            decay=0.36,
        )
        channel_rank *= np.maximum(rank, 0.42)
        channel_rank = np.clip(channel_rank, 0.0, 1.0)

    terrace = (
        rank
        * depth_gate
        * (2.0 + 8.5 * rank)
        * (0.76 + 0.24 * ridged)
        * (0.82 + 0.18 * passive_margin.astype(np.float64))
    )
    canyon = (
        -channel_rank
        * depth_gate
        * (4.0 + 15.0 * np.maximum(-texture, 0.0))
        * (0.55 + 0.30 * active_margin.astype(np.float64))
    )
    delta = np.clip(terrace + canyon, -22.0, 24.0)
    delta *= _polar_detail_damping(grid)
    delta[~domain] = 0.0
    return delta


def _shelf_slope_axis_rank(
    grid: SphereGrid,
    domain: np.ndarray,
    score: np.ndarray,
    shelf_break_rank: np.ndarray,
) -> np.ndarray:
    domain = np.asarray(domain, dtype=bool)
    score = np.asarray(score, dtype=np.float64)
    shelf_break = np.asarray(shelf_break_rank, dtype=np.float64)
    rank = np.zeros(score.shape, dtype=np.float64)
    values = score[domain]
    if not values.size:
        return rank
    cut = max(0.46, float(np.nanpercentile(values, 78.0)))
    support_cut = max(0.34, cut * 0.64)
    seed = domain & (score >= cut)
    support = (
        domain
        & (score >= support_cut)
        & (shelf_break >= 0.48)
    )
    if not np.any(seed & support):
        return rank
    seed &= support
    bridged = _bridge_sparse_line_mask(grid, seed, domain=support, passes=2)
    neighbor_fraction = _neighbor_mask_fraction(grid, bridged)
    local_axis = support & (
        (score >= max(0.50, cut * 0.84))
        | (neighbor_fraction >= 0.25)
    )
    axis = _bridge_sparse_line_mask(
        grid,
        bridged | local_axis,
        domain=support,
        passes=1,
    )
    rank = _ranked_halo(
        grid,
        axis,
        domain=support,
        passes=1,
        decay=0.42,
    )
    rank *= np.maximum(score, 0.46)
    rank[rank < 0.10] = 0.0
    return np.clip(rank, 0.0, 1.0)


def _marine_deep_ocean_fabric_delta(
    grid: SphereGrid,
    rel: np.ndarray,
    refined_rel: np.ndarray,
    ocean_depth: np.ndarray,
    crust_age: np.ndarray,
    ocean_distance: np.ndarray,
    fracture_zone: np.ndarray,
    transform: np.ndarray,
    abyssal_plain: np.ndarray,
    ridge: np.ndarray,
    trench: np.ndarray,
    process_island_promoted: np.ndarray,
    low_polar: np.ndarray,
    detail_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Add low-amplitude deep-ocean fabric without making new shoals."""
    rel = np.asarray(rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    ocean_depth = np.asarray(ocean_depth)
    crust_age = np.asarray(crust_age, dtype=np.float64)
    d = np.asarray(ocean_distance, dtype=np.float64)
    fracture_zone = np.asarray(fracture_zone, dtype=bool)
    transform = np.asarray(transform, dtype=bool)
    abyssal_plain = np.asarray(abyssal_plain, dtype=bool)
    ridge = np.asarray(ridge, dtype=bool)
    trench = np.asarray(trench, dtype=bool)
    promoted = np.asarray(process_island_promoted, dtype=bool)
    low_polar = np.asarray(low_polar, dtype=bool)
    ocean = rel < 0.0
    deep_domain = (
        ocean
        & low_polar
        & ~ridge
        & ~trench
        & ~promoted
        & (d >= 5.0)
        & (rel < -2400.0)
        & (rel > -5850.0)
        & (refined_rel < -1200.0)
    )
    delta = np.zeros(rel.shape, dtype=np.float64)
    fracture_rank = np.zeros(rel.shape, dtype=np.float64)
    plain_rank = np.zeros(rel.shape, dtype=np.float64)
    if not np.any(deep_domain):
        return delta, fracture_rank, plain_rank

    old_ocean = np.clip(crust_age / 185.0, 0.0, 1.0)
    depth_gate = np.clip((-refined_rel - 2200.0) / 2600.0, 0.0, 1.0)
    texture = _signed_texture(grid, int(detail_seed) + 4253, scale=6)
    ridged = _ridged_texture(grid, int(detail_seed) + 4271)

    fracture_seed = deep_domain & (fracture_zone | transform)
    fracture_exclusion = np.zeros(rel.shape, dtype=bool)
    if np.any(fracture_seed):
        neighbor_count = np.zeros(rel.shape, dtype=np.int16)
        for cell in np.flatnonzero(fracture_seed):
            nbs = grid.neighbors[int(cell)]
            neighbor_count[int(cell)] = int(np.count_nonzero(fracture_seed[nbs]))
        line_candidate = fracture_seed & (neighbor_count <= 3)
        if np.count_nonzero(line_candidate) < max(8, int(0.15 * np.count_nonzero(fracture_seed))):
            line_candidate = fracture_seed & (neighbor_count <= 4)
        thin_seed = line_candidate if np.any(line_candidate) else fracture_seed
        thin_seed = _bridge_sparse_line_mask(
            grid,
            thin_seed,
            domain=deep_domain,
            passes=1,
        )
        fracture_exclusion = _dilate_mask(grid, thin_seed, passes=1) & deep_domain
        fracture_rank[thin_seed] = np.clip(
            0.70 + 0.24 * old_ocean + 0.18 * depth_gate,
            0.0,
            1.0,
        )[thin_seed]

    plain_domain = (
        deep_domain
        & abyssal_plain
        & (ocean_depth == OCEAN_DEPTH_ABYSS)
        & (d >= 7.0)
        & ~fracture_exclusion
    )
    if np.any(plain_domain):
        plain_score = np.clip(
            (0.48 + 0.28 * old_ocean + 0.24 * depth_gate)
            * (0.78 + 0.22 * ridged),
            0.0,
            1.0,
        )
        values = plain_score[plain_domain]
        cut = max(0.62, float(np.nanpercentile(values, 72.0)))
        plain_seed = plain_domain & (plain_score >= cut)
        if np.any(plain_seed):
            plain_rank = _ranked_halo(
                grid,
                plain_seed,
                domain=plain_domain,
                passes=1,
                decay=0.28,
            )
            plain_rank *= plain_score

    if np.any(plain_rank > 0.0):
        local_mean = _local_domain_neighbor_mean(grid, refined_rel, deep_domain)
        flatten = np.clip(local_mean - refined_rel, -16.0, 16.0)
        delta += 0.26 * plain_rank * flatten
        delta += plain_rank * depth_gate * (1.6 * texture)

    if np.any(fracture_rank > 0.0):
        core = np.clip(fracture_rank, 0.0, 1.0)
        delta -= core * depth_gate * (8.0 + 11.0 * (0.55 + 0.45 * ridged))

    delta = np.clip(delta, -24.0, 12.0)
    delta *= _polar_detail_damping(grid)
    delta[~deep_domain] = 0.0
    return (
        delta,
        np.clip(fracture_rank, 0.0, 1.0),
        np.clip(plain_rank, 0.0, 1.0),
    )


def _marine_island_atoll_microrelief_delta(
    grid: SphereGrid,
    rel: np.ndarray,
    refined_rel: np.ndarray,
    islet_microshape_rank: np.ndarray,
    atoll_microshape_rank: np.ndarray,
    reef_rim_rank: np.ndarray,
    atoll_lagoon_rank: np.ndarray,
    fringing_reef_rank: np.ndarray,
    process_island_promoted: np.ndarray,
    atoll_islet_promoted: np.ndarray,
    ocean: np.ndarray,
    ridge_texture: np.ndarray,
    low_polar: np.ndarray,
    detail_seed: int,
) -> np.ndarray:
    """Add small terrain relief for island/atoll microshapes.

    Non-promoted ocean candidates are allowed to become shallower, but the
    caller's final clamp keeps them ocean.  Promoted cells keep their explicit
    opt-in land contract.
    """
    rel = np.asarray(rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    ocean = np.asarray(ocean, dtype=bool)
    low_polar = np.asarray(low_polar, dtype=bool)
    promoted = np.asarray(process_island_promoted, dtype=bool)
    atoll_promoted = np.asarray(atoll_islet_promoted, dtype=bool)
    islet = np.sqrt(np.clip(islet_microshape_rank, 0.0, 1.0))
    atoll = np.sqrt(np.clip(atoll_microshape_rank, 0.0, 1.0))
    rim = np.sqrt(np.clip(reef_rim_rank, 0.0, 1.0))
    lagoon = np.sqrt(np.clip(atoll_lagoon_rank, 0.0, 1.0))
    fringing = np.sqrt(np.clip(fringing_reef_rank, 0.0, 1.0))
    ridge = np.clip(np.asarray(ridge_texture, dtype=np.float64), 0.0, 1.0)
    texture = _signed_texture(grid, int(detail_seed) + 3119, scale=5)
    shallow_gate = np.clip((1250.0 + refined_rel) / 1250.0, 0.0, 1.0)
    seed = (
        (islet > 0.0)
        | (atoll > 0.0)
        | promoted
        | atoll_promoted
    )
    halo_seed = _dilate_mask(grid, seed, passes=1) if np.any(seed) else seed
    domain = (
        ocean
        & low_polar
        & (shallow_gate > 0.08)
        & (
            (islet > 0.0)
            | (atoll > 0.0)
            | (rim > 0.0)
            | (lagoon > 0.0)
            | (fringing > 0.0)
            | promoted
            | halo_seed
        )
    )
    delta = np.zeros(rel.shape, dtype=np.float64)
    if not np.any(domain):
        return delta
    core_seed = (
        (islet >= 0.70)
        | (atoll >= 0.60)
        | promoted
        | atoll_promoted
    )
    soft_halo = (
        _ranked_halo(grid, core_seed, domain=domain, passes=1, decay=0.24)
        if np.any(core_seed)
        else np.zeros(grid.n, dtype=np.float64)
    )
    islet_axis = np.maximum(islet, 0.30 * soft_halo)
    atoll_axis = np.maximum(
        atoll,
        0.28 * soft_halo * np.maximum(atoll, atoll_promoted.astype(np.float64)),
    )
    reef_axis = np.maximum(islet_axis, atoll_axis)

    islet_buildup = (
        islet_axis
        * shallow_gate
        * (3.0 + 8.0 * ridge + 2.5 * np.maximum(texture, 0.0))
    )
    islet_buildup[promoted] += (
        (2.0 + 5.5 * np.maximum(islet_axis[promoted], 0.55))
        if np.any(promoted)
        else 0.0
    )
    reef_buildup = (
        np.maximum(rim, 0.72 * fringing)
        * np.maximum(reef_axis, 0.18 * soft_halo)
        * shallow_gate
        * (1.6 + 4.5 * (0.58 + 0.42 * ridge))
    )
    reef_buildup[(reef_axis <= 0.0) & (soft_halo <= 0.0)] = 0.0
    atoll_rim_buildup = (
        atoll_axis
        * np.maximum(rim, 0.38)
        * (1.0 - 0.46 * np.clip(lagoon, 0.0, 1.0))
        * shallow_gate
        * (2.4 + 6.2 * (0.55 + 0.45 * ridge))
    )
    lagoon_support = np.maximum(lagoon * np.maximum(atoll_axis, 0.35), 0.32 * atoll_axis)
    lagoon_cut = (
        -lagoon_support
        * shallow_gate
        * (4.0 + 10.0 * (0.52 + 0.48 * ridge))
    )
    lagoon_cut[atoll_promoted] *= 0.70

    delta[domain] = (
        islet_buildup[domain]
        + reef_buildup[domain]
        + atoll_rim_buildup[domain]
        + lagoon_cut[domain]
    )
    delta = np.clip(delta, -18.0, 20.0)
    delta *= _polar_detail_damping(grid)
    delta[~domain] = 0.0
    return delta


def _marine_submarine_highland_morphology(
    grid: SphereGrid,
    rel: np.ndarray,
    potential_rel: np.ndarray,
    seamount: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    abyssal_hill: np.ndarray,
    ridge: np.ndarray,
    trench: np.ndarray,
    ocean_depth: np.ndarray,
    ocean_distance: np.ndarray,
    ridge_texture: np.ndarray,
    fine_texture: np.ndarray,
    low_polar: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Expose small process-backed submarine highs without changing topology."""
    rel = np.asarray(rel, dtype=np.float64)
    potential_rel = np.asarray(potential_rel, dtype=np.float64)
    ocean = rel < 0.0
    low_polar = np.asarray(low_polar, dtype=bool)
    trench = np.asarray(trench, dtype=bool)
    ridge = np.asarray(ridge, dtype=bool)
    seamount = np.asarray(seamount, dtype=bool)
    plateau = np.asarray(plateau, dtype=bool)
    microcontinent = np.asarray(microcontinent, dtype=bool)
    abyssal_hill = np.asarray(abyssal_hill, dtype=bool)
    ridge_texture = np.asarray(ridge_texture, dtype=np.float64)
    fine_texture = np.asarray(fine_texture, dtype=np.float64)
    ocean_depth = np.asarray(ocean_depth)
    ocean_distance = np.asarray(ocean_distance, dtype=np.float64)

    seamount_peak_rank = np.zeros(grid.n, dtype=np.float64)
    seamount_apron_rank = np.zeros(grid.n, dtype=np.float64)
    plateau_edge_rank = np.zeros(grid.n, dtype=np.float64)
    abyssal_hill_rank = np.zeros(grid.n, dtype=np.float64)

    shallow = np.clip((3200.0 + rel) / 3200.0, 0.0, 1.0)
    uplift = np.clip((potential_rel - rel) / 620.0, 0.0, 1.0)
    seamount_domain = ocean & seamount & low_polar & ~trench & (rel > -3400.0)
    if seamount_domain.any():
        score = np.clip(
            0.36 * shallow
            + 0.34 * ridge_texture
            + 0.22 * uplift
            + 0.08 * np.clip((72.0 - np.abs(grid.lat)) / 72.0, 0.0, 1.0),
            0.0,
            1.0,
        )
        seamount_peak_rank, chain_axis = _seamount_chain_peak_and_axis_rank(
            grid,
            seamount_domain,
            score,
        )
        apron_domain = (
            ocean
            & low_polar
            & ~trench
            & ~ridge
            & (rel > -4300.0)
            & (ocean_distance >= 0.0)
        )
        apron_seed = (chain_axis > 0.0) | (seamount_peak_rank > 0.0)
        seamount_apron_rank = _ranked_halo(
            grid,
            apron_seed,
            domain=apron_domain,
            passes=2,
            decay=0.44,
        )
        seamount_apron_rank *= np.maximum(
            0.42 * chain_axis,
            np.where(seamount_peak_rank > 0.0, 0.72, 0.54),
        )
        seamount_apron_rank[seamount_peak_rank > 0.0] = np.maximum(
            seamount_apron_rank[seamount_peak_rank > 0.0],
            0.48 * seamount_peak_rank[seamount_peak_rank > 0.0],
        )

    plateau_domain = ocean & (plateau | microcontinent) & low_polar & ~trench
    if plateau_domain.any():
        edge = np.zeros(grid.n, dtype=bool)
        for cell in np.flatnonzero(plateau_domain):
            nbs = grid.neighbors[int(cell)]
            if nbs.size and np.any(~plateau_domain[nbs] & ocean[nbs]):
                edge[int(cell)] = True
        if edge.any():
            edge_domain = (
                ocean
                & low_polar
                & ~trench
                & (rel > -3600.0)
                & (ocean_distance >= 0.0)
            )
            plateau_edge_rank = _plateau_escarpment_segment_rank(
                grid,
                edge,
                edge_domain,
                plateau,
                microcontinent,
                ridge_texture,
                fine_texture,
            )
            plateau_edge_rank[plateau_edge_rank < 0.12] = 0.0

    abyssal_domain = (
        ocean
        & abyssal_hill
        & low_polar
        & ~trench
        & ~ridge
        & (ocean_depth == OCEAN_DEPTH_ABYSS)
        & (rel < -2200.0)
        & (rel > -5600.0)
    )
    if abyssal_domain.any():
        texture_strength = np.clip(0.52 * ridge_texture + 0.48 * np.abs(fine_texture),
                                   0.0,
                                   1.0)
        cut = max(0.74, float(np.nanpercentile(texture_strength[abyssal_domain], 94.0)))
        abyssal_seed = abyssal_domain & (texture_strength >= cut)
        if abyssal_seed.any():
            abyssal_hill_rank = _ranked_halo(
                grid,
                abyssal_seed,
                domain=abyssal_domain,
                passes=1,
                decay=0.34,
            )
            abyssal_hill_rank *= np.maximum(texture_strength, 0.52)
            abyssal_hill_rank[abyssal_hill_rank < 0.20] = 0.0

    delta = np.zeros(grid.n, dtype=np.float64)
    peak_cells = seamount_peak_rank > 0.0
    peak_target = -180.0 - 520.0 * (1.0 - seamount_peak_rank)
    delta[peak_cells] += np.minimum(
        135.0 * seamount_peak_rank[peak_cells],
        np.maximum(0.0, peak_target[peak_cells] - rel[peak_cells]) * 0.24,
    )
    apron_cells = seamount_apron_rank > 0.05
    delta[apron_cells] += 24.0 * seamount_apron_rank[apron_cells]

    plateau_cells = plateau_edge_rank > 0.05
    delta[plateau_cells] -= 22.0 * plateau_edge_rank[plateau_cells]

    abyssal_cells = abyssal_hill_rank > 0.0
    delta[abyssal_cells] += 18.0 * abyssal_hill_rank[abyssal_cells]
    delta[trench | ridge | ~ocean] = 0.0
    return (
        delta,
        np.clip(seamount_peak_rank, 0.0, 1.0),
        np.clip(seamount_apron_rank, 0.0, 1.0),
        np.clip(plateau_edge_rank, 0.0, 1.0),
        np.clip(abyssal_hill_rank, 0.0, 1.0),
    )


def _seamount_chain_peak_and_axis_rank(
    grid: SphereGrid,
    domain: np.ndarray,
    score: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    domain = np.asarray(domain, dtype=bool)
    score = np.asarray(score, dtype=np.float64)
    peak_rank = np.zeros(score.shape, dtype=np.float64)
    axis_rank = np.zeros(score.shape, dtype=np.float64)
    values = score[domain]
    if not values.size:
        return peak_rank, axis_rank

    core_cut = max(0.54, float(np.nanpercentile(values, 72.0)))
    support_cut = max(0.40, float(np.nanpercentile(values, 45.0)))
    support = domain & (score >= support_cut)
    candidates = np.flatnonzero(domain & (score >= core_cut))
    if candidates.size == 0:
        return peak_rank, axis_rank

    max_peaks = min(
        int(candidates.size),
        max(8, int(round(0.0011 * grid.n))),
    )
    peak_seed = _spaced_rank_seed(
        grid,
        candidates,
        score,
        max_count=max_peaks,
        spacing_passes=1,
    )
    if not np.any(peak_seed):
        return peak_rank, axis_rank

    axis_seed = _bridge_sparse_line_mask(
        grid,
        peak_seed,
        domain=support,
        passes=2,
    )
    axis_seed |= peak_seed
    axis_rank = _ranked_halo(
        grid,
        axis_seed,
        domain=support,
        passes=1,
        decay=0.38,
    )
    axis_rank *= np.maximum(score, 0.48)
    axis_rank[axis_rank < 0.10] = 0.0

    peak_rank[peak_seed] = score[peak_seed]
    peak_rank[peak_seed] = np.maximum(peak_rank[peak_seed], score[peak_seed])
    return np.clip(peak_rank, 0.0, 1.0), np.clip(axis_rank, 0.0, 1.0)


def _plateau_escarpment_segment_rank(
    grid: SphereGrid,
    edge: np.ndarray,
    edge_domain: np.ndarray,
    plateau: np.ndarray,
    microcontinent: np.ndarray,
    ridge_texture: np.ndarray,
    fine_texture: np.ndarray,
) -> np.ndarray:
    edge = np.asarray(edge, dtype=bool)
    edge_domain = np.asarray(edge_domain, dtype=bool)
    plateau = np.asarray(plateau, dtype=bool)
    microcontinent = np.asarray(microcontinent, dtype=bool)
    ridge = np.clip(np.asarray(ridge_texture, dtype=np.float64), 0.0, 1.0)
    fine = np.asarray(fine_texture, dtype=np.float64)
    rank = np.zeros(ridge.shape, dtype=np.float64)
    if not np.any(edge):
        return rank

    texture_score = np.clip(
        0.58 * ridge
        + 0.26 * np.abs(fine)
        + 0.16 * _neighbor_mask_fraction(grid, edge),
        0.0,
        1.0,
    )
    values = texture_score[edge]
    if not values.size:
        return rank
    cut = float(np.nanpercentile(values, 58.0))
    segment_seed = edge & (texture_score >= cut)
    if not np.any(segment_seed):
        cells = np.flatnonzero(edge)
        segment_seed[int(cells[np.nanargmax(texture_score[cells])])] = True

    segment_seed = _bridge_sparse_line_mask(
        grid,
        segment_seed,
        domain=edge,
        passes=1,
    )
    edge_band = _dilate_mask(grid, edge, passes=1) & edge_domain
    rank = _ranked_halo(
        grid,
        segment_seed,
        domain=edge_band,
        passes=1,
        decay=0.40,
    )
    class_weight = np.where(plateau, 0.92, np.where(microcontinent, 0.78, 0.48))
    rank *= class_weight * (0.72 + 0.28 * texture_score)
    plateau_seed = segment_seed & plateau
    micro_seed = segment_seed & microcontinent
    rank[plateau_seed] = np.maximum(rank[plateau_seed], 0.82)
    rank[micro_seed] = np.maximum(rank[micro_seed], 0.72)
    rank[rank < 0.12] = 0.0
    return np.clip(rank, 0.0, 1.0)


def _spaced_rank_seed(
    grid: SphereGrid,
    candidates: np.ndarray,
    values: np.ndarray,
    *,
    max_count: int,
    spacing_passes: int,
) -> np.ndarray:
    vals = np.asarray(values, dtype=np.float64)
    cells = np.asarray(candidates, dtype=np.int64)
    seed = np.zeros(vals.shape, dtype=bool)
    if cells.size == 0 or int(max_count) <= 0:
        return seed
    ordered = cells[np.argsort(vals[cells])[::-1]]
    blocked = np.zeros(vals.shape, dtype=bool)
    selected = 0
    for cell in ordered:
        cell = int(cell)
        if blocked[cell]:
            continue
        seed[cell] = True
        selected += 1
        if selected >= int(max_count):
            break
        if int(spacing_passes) > 0:
            _block_local_neighborhood(
                grid,
                blocked,
                cell,
                passes=int(spacing_passes),
            )
    return seed


def _dilate_mask(grid: SphereGrid, mask: np.ndarray, *, passes: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(0, int(passes))):
        expanded = out.copy()
        for cell in np.flatnonzero(out):
            expanded[grid.neighbors[int(cell)]] = True
        out = expanded
    return out


def _mask_components(grid: SphereGrid, mask: np.ndarray) -> list[np.ndarray]:
    domain = np.asarray(mask, dtype=bool)
    visited = np.zeros(grid.n, dtype=bool)
    out: list[np.ndarray] = []
    for start in np.flatnonzero(domain):
        start = int(start)
        if visited[start]:
            continue
        cells: list[int] = []
        q: deque[int] = deque([start])
        visited[start] = True
        while q:
            cell = int(q.popleft())
            cells.append(cell)
            for nb_raw in grid.neighbors[cell]:
                nb = int(nb_raw)
                if domain[nb] and not visited[nb]:
                    visited[nb] = True
                    q.append(nb)
        out.append(np.asarray(cells, dtype=np.int64))
    return out


def _neighbor_mask_count(grid: SphereGrid, mask: np.ndarray) -> np.ndarray:
    domain = np.asarray(mask, dtype=bool)
    count = np.zeros(grid.n, dtype=np.int16)
    for cell in range(grid.n):
        nbs = grid.neighbors[int(cell)]
        if nbs.size:
            count[int(cell)] = int(np.count_nonzero(domain[nbs]))
    return count


def _unsupported_inland_shallow_seaway_tuning(
    grid: SphereGrid,
    rel: np.ndarray,
    ocean: np.ndarray,
    ocean_distance: np.ndarray,
    ocean_depth: np.ndarray,
    shelf_width: np.ndarray,
    gateway_id: np.ndarray,
    gateway_system_id: np.ndarray,
    process_support: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Narrow unsupported shallow gateway belts into lowland plus a small thalweg."""
    rel = np.asarray(rel, dtype=np.float64)
    ocean = np.asarray(ocean, dtype=bool)
    d = np.asarray(ocean_distance, dtype=np.float64)
    depth = np.asarray(ocean_depth, dtype=np.int32)
    shelf = np.asarray(shelf_width, dtype=np.float64)
    gateway = (
        (np.asarray(gateway_id, dtype=np.int32) >= 0)
        | (np.asarray(gateway_system_id, dtype=np.int32) >= 0)
        | (depth == int(OCEAN_DEPTH_RESTRICTED))
    )
    support = np.asarray(process_support, dtype=bool)
    delta = np.zeros(grid.n, dtype=np.float64)
    tuned = np.zeros(grid.n, dtype=bool)
    rank = np.zeros(grid.n, dtype=np.float64)
    landback = np.zeros(grid.n, dtype=bool)
    if not np.any(ocean & gateway):
        return delta, tuned, rank, landback

    land = ~ocean
    land_neighbor_count = _neighbor_mask_count(grid, land)
    ocean_neighbor_count = _neighbor_mask_count(grid, ocean)
    support_halo = _dilate_mask(grid, support, passes=1)
    landlocked_shallow_fallback = (
        (land_neighbor_count >= 3)
        & (ocean_neighbor_count >= 2)
        & (d <= 2.75)
        & (shelf <= 4.0)
    )
    constricted = (
        ocean
        & (gateway | landlocked_shallow_fallback)
        & np.isfinite(d)
        & (d <= 3.35)
        & (rel > -980.0)
        & (rel < -8.0)
        & (land_neighbor_count >= 2)
        & (ocean_neighbor_count >= 2)
        & ~support_halo
    )
    if not np.any(constricted):
        return delta, tuned, rank, landback

    min_cells = max(4, int(round(0.000045 * grid.n)))
    for comp in _mask_components(grid, constricted):
        if comp.size < min_cells:
            continue
        land_lock = float(np.mean(land_neighbor_count[comp] >= 2))
        gateway_fraction = float(np.mean(gateway[comp]))
        support_fraction = float(np.mean(support_halo[comp]))
        median_rel = float(np.nanmedian(rel[comp]))
        fallback_fraction = float(np.mean(landlocked_shallow_fallback[comp]))
        if (
            land_lock < 0.62
            or (gateway_fraction < 0.35 and fallback_fraction < 0.72)
            or support_fraction > 0.08
            or median_rel <= -980.0
        ):
            continue
        local_center = np.clip(
            (ocean_neighbor_count[comp].astype(np.float64)
             - land_neighbor_count[comp].astype(np.float64) + 3.0) / 6.0,
            0.0,
            1.0,
        )
        pass_factor = np.clip(d[comp] / 3.35, 0.0, 1.0)
        restricted_bonus = (depth[comp] == int(OCEAN_DEPTH_RESTRICTED)).astype(np.float64)
        shelf_factor = np.clip(shelf[comp] / 4.0, 0.0, 1.0)
        comp_rank = np.clip(
            0.22
            + 0.32 * local_center
            + 0.28 * pass_factor
            + 0.12 * restricted_bonus
            + 0.06 * shelf_factor,
            0.0,
            1.0,
        )
        keep_fraction = np.clip(
            0.05 + 0.10 * np.mean(local_center) + 0.04 * (1.0 - land_lock),
            0.04,
            0.18,
        )
        keep_count = max(1, min(comp.size, int(round(float(comp.size) * keep_fraction))))
        order = np.argsort(comp_rank)
        thalweg = np.zeros(comp.size, dtype=bool)
        thalweg[order[-keep_count:]] = True
        thalweg_cells = comp[thalweg]
        landback_cells = comp[~thalweg]

        target = -(
            420.0
            + 360.0 * comp_rank
            + 140.0 * pass_factor
            + 160.0 * restricted_bonus
        )
        deepen = np.minimum(0.0, 0.80 * (target[thalweg] - rel[thalweg_cells]))
        strong = deepen < -12.0
        if np.any(strong):
            cells = thalweg_cells[strong]
            delta[cells] = deepen[strong]
            tuned[cells] = True
            rank[cells] = comp_rank[thalweg][strong]

        if landback_cells.size:
            relief_texture = _signed_texture(grid, 61001 + int(comp.size), scale=3)
            lowland_target = (
                12.0
                + 38.0 * np.clip(1.0 - comp_rank[~thalweg], 0.0, 1.0)
                + 10.0 * relief_texture[landback_cells]
            )
            lift = np.maximum(0.0, lowland_target - rel[landback_cells])
            close = lift > 10.0
            if np.any(close):
                cells = landback_cells[close]
                delta[cells] = lift[close]
                tuned[cells] = True
                landback[cells] = True
                rank[cells] = np.clip(1.0 - comp_rank[~thalweg][close], 0.0, 1.0)

    return delta, tuned, rank, landback


def _final_unsupported_shallow_corridor_landback(
    grid: SphereGrid,
    parent_rel: np.ndarray,
    refined_rel: np.ndarray,
    parent_ocean: np.ndarray,
    process_support: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Close residual unsupported shallow water trapped between land cells."""
    parent_rel = np.asarray(parent_rel, dtype=np.float64)
    refined_rel = np.asarray(refined_rel, dtype=np.float64)
    parent_ocean = np.asarray(parent_ocean, dtype=bool)
    support = _dilate_mask(grid, np.asarray(process_support, dtype=bool), passes=1)
    current_land = refined_rel >= 0.0
    current_ocean = ~current_land
    land_neighbor_count = _neighbor_mask_count(grid, current_land)
    ocean_neighbor_count = _neighbor_mask_count(grid, current_ocean)
    candidate = (
        parent_ocean
        & current_ocean
        & (refined_rel > -820.0)
        & (parent_rel > -1100.0)
        & (land_neighbor_count >= 3)
        & (ocean_neighbor_count >= 2)
        & ~support
    )
    delta = np.zeros(grid.n, dtype=np.float64)
    landback = np.zeros(grid.n, dtype=bool)
    if not np.any(candidate):
        return delta, landback

    for comp in _mask_components(grid, candidate):
        if comp.size < 2:
            continue
        land_lock = float(np.mean(land_neighbor_count[comp] >= 3))
        if land_lock < 0.72:
            continue
        local_relief = _signed_texture(grid, 71003 + int(comp.size), scale=4)
        target = 10.0 + 32.0 * land_lock + 8.0 * local_relief[comp]
        lift = np.maximum(0.0, target - refined_rel[comp])
        close = lift > 8.0
        if not np.any(close):
            continue
        cells = comp[close]
        delta[cells] = lift[close]
        landback[cells] = True
    return delta, landback


def _bridge_sparse_line_mask(
    grid: SphereGrid,
    seed: np.ndarray,
    *,
    domain: np.ndarray,
    passes: int,
) -> np.ndarray:
    line = np.asarray(seed, dtype=bool).copy() & np.asarray(domain, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    for _ in range(max(0, int(passes))):
        search = _dilate_mask(grid, line, passes=1) & domain & ~line
        if not np.any(search):
            break
        add = np.zeros(grid.n, dtype=bool)
        for cell in np.flatnonzero(search):
            nbs = grid.neighbors[int(cell)]
            if np.count_nonzero(line[nbs]) == 2:
                add[int(cell)] = True
        if not np.any(add):
            break
        line |= add
    return line


def _neighbor_mask_fraction(grid: SphereGrid, mask: np.ndarray) -> np.ndarray:
    domain = np.asarray(mask, dtype=bool)
    fraction = np.zeros(grid.n, dtype=np.float64)
    for cell in range(grid.n):
        nbs = grid.neighbors[int(cell)]
        if nbs.size:
            fraction[int(cell)] = float(np.mean(domain[nbs]))
    return fraction


def _local_domain_neighbor_mean(
    grid: SphereGrid,
    values: np.ndarray,
    domain: np.ndarray,
) -> np.ndarray:
    vals = np.asarray(values, dtype=np.float64)
    domain = np.asarray(domain, dtype=bool)
    out = vals.copy()
    for cell in np.flatnonzero(domain):
        nbs = grid.neighbors[int(cell)]
        local = nbs[domain[nbs]]
        if local.size:
            out[int(cell)] = float(np.nanmean(vals[local]))
    return out


def _ranked_halo(
    grid: SphereGrid,
    seed: np.ndarray,
    *,
    domain: np.ndarray,
    passes: int,
    decay: float,
) -> np.ndarray:
    rank = np.zeros(grid.n, dtype=np.float64)
    current = np.asarray(seed, dtype=bool) & np.asarray(domain, dtype=bool)
    if not current.any():
        return rank
    rank[current] = 1.0
    visited = current.copy()
    strength = 1.0
    for _ in range(max(0, int(passes))):
        strength *= float(decay)
        grown = np.zeros(grid.n, dtype=bool)
        for cell in np.flatnonzero(current):
            grown[grid.neighbors[int(cell)]] = True
        grown &= np.asarray(domain, dtype=bool) & ~visited
        if not grown.any():
            break
        rank[grown] = np.maximum(rank[grown], strength)
        visited |= grown
        current = grown
    return np.clip(rank, 0.0, 1.0)


def _field_key(fields: dict[str, str], field_name: str) -> str:
    return str(fields.get(field_name) or "field__" + field_name.replace(".", "_"))


def _int_field(
    context: _RefinementContext,
    field_name: str,
    *,
    default: int | np.ndarray,
) -> np.ndarray:
    key = _field_key(context.fields, field_name)
    if key not in context.arrays:
        if isinstance(default, np.ndarray):
            return default.astype(np.int32)
        return np.full(context.target_grid.n, int(default), dtype=np.int32)
    return np.asarray(context.arrays[key], dtype=np.float64).astype(np.int32)


def _float_field(
    context: _RefinementContext,
    field_name: str,
    *,
    default: float,
) -> np.ndarray:
    key = _field_key(context.fields, field_name)
    if key not in context.arrays:
        return np.full(context.target_grid.n, float(default), dtype=np.float64)
    return np.asarray(context.arrays[key], dtype=np.float64)


def _mask(context: _RefinementContext, key: str) -> np.ndarray:
    raw = context.arrays.get(key)
    if raw is None:
        return np.zeros(context.target_grid.n, dtype=bool)
    return np.asarray(raw, dtype=bool)


def _signed_texture(grid: SphereGrid, seed: int, *, scale: int) -> np.ndarray:
    phase = float(seed % 1000003) * 0.000173
    s = float(max(int(scale), 1))
    directions = np.asarray([
        [0.71, 0.17, 0.68],
        [-0.28, 0.83, 0.48],
        [0.58, -0.63, 0.52],
        [-0.74, -0.39, 0.55],
        [0.18, 0.94, -0.29],
        [0.91, -0.12, -0.40],
    ], dtype=np.float64)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    freqs = np.asarray([13.0, 19.0, 29.0, 41.0, 53.0, 67.0],
                       dtype=np.float64) * s
    amps = np.asarray([0.30, 0.23, 0.18, 0.14, 0.10, 0.05],
                      dtype=np.float64)
    dots = grid.xyz @ directions.T
    value = np.zeros(grid.n, dtype=np.float64)
    for i in range(directions.shape[0]):
        value += amps[i] * np.sin(freqs[i] * dots[:, i] + phase * (i + 1.0))
    return np.clip(value / np.sum(amps), -1.0, 1.0)


def _ridged_texture(grid: SphereGrid, seed: int) -> np.ndarray:
    tex = _signed_texture(grid, seed, scale=3)
    return np.clip(1.0 - np.abs(tex), 0.0, 1.0)


def _polar_detail_damping(grid: SphereGrid) -> np.ndarray:
    lat_abs = np.abs(np.asarray(grid.lat, dtype=np.float64))
    return np.clip((86.0 - lat_abs) / 12.0, 0.25, 1.0)


def _refined_metrics(
    config: SelectedSnapshotRefinementConfig,
    source_metrics_path: Path,
    source_arrays_path: Path,
    refined_arrays_path: Path,
    parent_metrics: dict[str, Any],
    parent_rel: np.ndarray,
    refined_rel: np.ndarray,
    detail_delta: np.ndarray,
    amplitude: np.ndarray,
    marine: _MarineMicrogeomorphology,
    micro: _Microgeomorphology,
    coastal: _CoastalMorphology,
    context: _RefinementContext,
) -> dict[str, Any]:
    parent_land = parent_rel >= 0.0
    refined_land = refined_rel >= 0.0
    land = parent_land
    ocean = ~land

    def pct(values: np.ndarray, q: float) -> float:
        values = np.asarray(values, dtype=np.float64)
        if values.size == 0:
            return 0.0
        return float(np.nanpercentile(values, q))

    def frac(mask: np.ndarray) -> float:
        return float(np.mean(np.asarray(mask, dtype=bool)))

    hydrology_nonzero = np.abs(micro.hydrology_delta[np.abs(micro.hydrology_delta) > 0.0])
    fluvial_microrelief_nonzero = np.abs(
        micro.fluvial_microrelief_delta[
            np.abs(micro.fluvial_microrelief_delta) > 0.0
        ]
    )
    lowland_alluvial_nonzero = np.abs(
        micro.lowland_alluvial_microrelief_delta[
            np.abs(micro.lowland_alluvial_microrelief_delta) > 0.0
        ]
    )
    marine_nonzero = np.abs(marine.marine_delta[np.abs(marine.marine_delta) > 0.0])
    shelf_slope_microrelief_nonzero = np.abs(
        marine.shelf_slope_microrelief_delta[
            np.abs(marine.shelf_slope_microrelief_delta) > 0.0
        ]
    )
    deep_ocean_fabric_nonzero = np.abs(
        marine.deep_ocean_fabric_delta[
            np.abs(marine.deep_ocean_fabric_delta) > 0.0
        ]
    )
    submarine_highland_nonzero = np.abs(
        marine.submarine_highland_delta[
            np.abs(marine.submarine_highland_delta) > 0.0
        ]
    )
    inland_seaway_tuning_nonzero = np.abs(
        marine.inland_seaway_tuning_delta[
            np.abs(marine.inland_seaway_tuning_delta) > 0.0
        ]
    )
    island_atoll_microrelief_nonzero = np.abs(
        marine.island_atoll_microrelief_delta[
            np.abs(marine.island_atoll_microrelief_delta) > 0.0
        ]
    )
    island_promotion_nonzero = np.abs(
        marine.process_island_promotion_delta[
            np.abs(marine.process_island_promotion_delta) > 0.0
        ]
    )
    coastal_nonzero = np.abs(coastal.coastal_delta[np.abs(coastal.coastal_delta) > 0.0])
    coastal_process_nonzero = np.abs(
        coastal.coastal_process_microrelief_delta[
            np.abs(coastal.coastal_process_microrelief_delta) > 0.0
        ]
    )
    coastal_depositional_nonzero = np.abs(
        coastal.coastal_depositional_microrelief_delta[
            np.abs(coastal.coastal_depositional_microrelief_delta) > 0.0
        ]
    )
    process_masks = {
        "orogen_hierarchy": _int_field(
            context, "terrain.orogenic_parent_hierarchy", default=0) > 0,
        "orogen_spine": _int_field(
            context, "terrain.orogenic_hierarchy_spine", default=0) > 0,
        "ridge": _mask(context, "boundary__ridge"),
        "trench": _mask(context, "boundary__trench"),
        "seamount_chain": _mask(context, "object__seamount_chain"),
        "oceanic_plateau": _mask(context, "object__oceanic_plateau"),
        "microcontinent": _mask(context, "object__microcontinent"),
        "island_arc": _mask(context, "object__island_arc"),
    }
    return {
        "schema": SCHEMA,
        "source_metrics": str(source_metrics_path),
        "source_arrays": str(source_arrays_path),
        "refined_arrays": str(refined_arrays_path),
        "source_cells": int(parent_metrics.get("cells", parent_rel.size)),
        "target_cells": int(config.target_cells),
        "width": int(config.width),
        "height": int(config.height),
        "sea_level_m": float(parent_metrics.get("sea_level_m", 0.0)),
        "detail_seed": int(config.detail_seed),
        "detail_strength": float(config.detail_strength),
        "allow_process_islands": bool(config.allow_process_islands),
        "render_groups": sorted(_normalize_render_groups(config.render_groups)),
        "land_fraction_parent": frac(parent_land),
        "land_fraction_refined": frac(refined_land),
        "land_ocean_sign_flip_fraction": frac(parent_land != refined_land),
        "detail_delta_abs_p50_m": pct(np.abs(detail_delta), 50.0),
        "detail_delta_abs_p95_m": pct(np.abs(detail_delta), 95.0),
        "detail_delta_land_abs_p95_m": pct(np.abs(detail_delta[land]), 95.0),
        "detail_delta_ocean_abs_p95_m": pct(np.abs(detail_delta[ocean]), 95.0),
        "detail_amplitude_p95_m": pct(amplitude, 95.0),
        "marine_delta_abs_p95_m": pct(np.abs(marine.marine_delta), 95.0),
        "marine_delta_nonzero_abs_p95_m": pct(marine_nonzero, 95.0),
        "shelf_slope_microrelief_delta_abs_p95_m": pct(
            np.abs(marine.shelf_slope_microrelief_delta), 95.0),
        "shelf_slope_microrelief_delta_nonzero_abs_p95_m": pct(
            shelf_slope_microrelief_nonzero, 95.0),
        "shelf_slope_microrelief_cell_fraction_ocean": (
            float(np.count_nonzero(
                (np.abs(marine.shelf_slope_microrelief_delta) > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "deep_ocean_fabric_delta_abs_p95_m": pct(
            np.abs(marine.deep_ocean_fabric_delta), 95.0),
        "deep_ocean_fabric_delta_nonzero_abs_p95_m": pct(
            deep_ocean_fabric_nonzero, 95.0),
        "deep_ocean_fabric_cell_fraction_ocean": (
            float(np.count_nonzero(
                (np.abs(marine.deep_ocean_fabric_delta) > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "submarine_highland_delta_nonzero_abs_p95_m": pct(
            submarine_highland_nonzero, 95.0),
        "inland_seaway_tuning_delta_nonzero_abs_p95_m": pct(
            inland_seaway_tuning_nonzero, 95.0),
        "inland_seaway_tuning_cell_fraction_ocean": (
            float(np.count_nonzero(marine.inland_seaway_tuning_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "inland_seaway_landback_cell_fraction_parent_ocean": (
            float(np.count_nonzero(marine.inland_seaway_landback_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "island_atoll_microrelief_delta_abs_p95_m": pct(
            np.abs(marine.island_atoll_microrelief_delta), 95.0),
        "island_atoll_microrelief_delta_nonzero_abs_p95_m": pct(
            island_atoll_microrelief_nonzero, 95.0),
        "island_atoll_microrelief_cell_fraction_ocean": (
            float(np.count_nonzero(
                (np.abs(marine.island_atoll_microrelief_delta) > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "process_island_promotion_delta_nonzero_abs_p95_m": pct(
            island_promotion_nonzero, 95.0),
        "hydrology_delta_abs_p95_m": pct(np.abs(micro.hydrology_delta), 95.0),
        "hydrology_delta_nonzero_abs_p95_m": pct(hydrology_nonzero, 95.0),
        "fluvial_microrelief_delta_abs_p95_m": pct(
            np.abs(micro.fluvial_microrelief_delta), 95.0),
        "fluvial_microrelief_delta_nonzero_abs_p95_m": pct(
            fluvial_microrelief_nonzero, 95.0),
        "fluvial_microrelief_cell_fraction_land": (
            float(np.count_nonzero((np.abs(micro.fluvial_microrelief_delta) > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "lowland_alluvial_microrelief_delta_abs_p95_m": pct(
            np.abs(micro.lowland_alluvial_microrelief_delta), 95.0),
        "lowland_alluvial_microrelief_delta_nonzero_abs_p95_m": pct(
            lowland_alluvial_nonzero, 95.0),
        "lowland_alluvial_microrelief_cell_fraction_land": (
            float(np.count_nonzero(
                (np.abs(micro.lowland_alluvial_microrelief_delta) > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "coastal_delta_abs_p95_m": pct(np.abs(coastal.coastal_delta), 95.0),
        "coastal_delta_nonzero_abs_p95_m": pct(coastal_nonzero, 95.0),
        "coastal_process_microrelief_delta_abs_p95_m": pct(
            np.abs(coastal.coastal_process_microrelief_delta), 95.0),
        "coastal_process_microrelief_delta_nonzero_abs_p95_m": pct(
            coastal_process_nonzero, 95.0),
        "coastal_process_microrelief_cell_fraction_ocean": (
            float(np.count_nonzero(
                (np.abs(coastal.coastal_process_microrelief_delta) > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "coastal_depositional_microrelief_delta_abs_p95_m": pct(
            np.abs(coastal.coastal_depositional_microrelief_delta), 95.0),
        "coastal_depositional_microrelief_delta_nonzero_abs_p95_m": pct(
            coastal_depositional_nonzero, 95.0),
        "coastal_depositional_microrelief_cell_fraction_land": (
            float(np.count_nonzero(
                (np.abs(coastal.coastal_depositional_microrelief_delta) > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "river_cell_fraction_land": (
            float(np.count_nonzero((micro.river_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "major_river_cell_fraction_land": (
            float(np.count_nonzero((micro.river_rank >= 0.55) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "river_path_cell_fraction_land": (
            float(np.count_nonzero((micro.river_path_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "drainage_basin_count": int(np.max(micro.drainage_basin_id)),
        "basin_trunk_cell_fraction_land": (
            float(np.count_nonzero((micro.basin_trunk_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "floodplain_cell_fraction_land": (
            float(np.count_nonzero((micro.floodplain_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "meander_belt_cell_fraction_land": (
            float(np.count_nonzero((micro.meander_belt_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "meander_scroll_cell_fraction_land": (
            float(np.count_nonzero((micro.meander_scroll_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "floodplain_swale_cell_fraction_land": (
            float(np.count_nonzero((micro.floodplain_swale_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "basin_trunk_distance_p95_passes": pct(
            micro.land_coast_distance[
                (micro.basin_trunk_rank > 0.0) & (micro.land_coast_distance >= 0.0)
            ],
            95.0,
        ),
        "lake_cell_fraction_land": (
            float(np.count_nonzero(micro.lake_mask & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "lake_basin_cell_fraction_land": (
            float(np.count_nonzero((micro.lake_basin_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "lake_shoreline_cell_fraction_land": (
            float(np.count_nonzero((micro.lake_shoreline_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "alluvial_fan_rank_cell_fraction_land": (
            float(np.count_nonzero((micro.alluvial_fan_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "lowland_plain_rank_cell_fraction_land": (
            float(np.count_nonzero((micro.lowland_plain_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "piedmont_apron_rank_cell_fraction_land": (
            float(np.count_nonzero((micro.piedmont_apron_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "delta_cell_fraction_ocean": (
            float(np.count_nonzero(micro.delta_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "delta_fan_cell_fraction_ocean": (
            float(np.count_nonzero((micro.delta_fan_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "delta_plain_cell_fraction_land": (
            float(np.count_nonzero(micro.delta_plain_mask & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "reef_atoll_cell_fraction_ocean": (
            float(np.count_nonzero(marine.reef_atoll_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "marine_shoal_cell_fraction_ocean": (
            float(np.count_nonzero(marine.marine_shoal_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "seamount_shoal_cell_fraction_ocean": (
            float(np.count_nonzero(marine.seamount_shoal_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "oceanic_plateau_shoal_cell_fraction_ocean": (
            float(np.count_nonzero(marine.oceanic_plateau_shoal_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "microcontinent_shoal_cell_fraction_ocean": (
            float(np.count_nonzero(marine.microcontinent_shoal_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "island_arc_shoal_cell_fraction_ocean": (
            float(np.count_nonzero(marine.island_arc_shoal_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "shelf_break_cell_fraction_ocean": (
            float(np.count_nonzero((marine.shelf_break_rank >= 0.68) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "process_island_candidate_cell_fraction_ocean": (
            float(np.count_nonzero(marine.process_island_candidate_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "island_candidate_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.island_candidate_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "atoll_candidate_cell_fraction_ocean": (
            float(np.count_nonzero(marine.atoll_candidate_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "atoll_candidate_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.atoll_candidate_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "reef_rim_cell_fraction_ocean": (
            float(np.count_nonzero((marine.reef_rim_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "atoll_lagoon_cell_fraction_ocean": (
            float(np.count_nonzero((marine.atoll_lagoon_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "fringing_reef_cell_fraction_ocean": (
            float(np.count_nonzero((marine.fringing_reef_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "process_island_promotion_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.process_island_promotion_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "islet_microshape_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.islet_microshape_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "atoll_microshape_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.atoll_microshape_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "process_island_promoted_cell_fraction_parent_ocean": (
            float(np.count_nonzero(marine.process_island_promoted_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "atoll_islet_promoted_cell_fraction_parent_ocean": (
            float(np.count_nonzero(marine.atoll_islet_promoted_mask & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "seamount_peak_cell_fraction_ocean": (
            float(np.count_nonzero((marine.seamount_peak_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "seamount_apron_cell_fraction_ocean": (
            float(np.count_nonzero((marine.seamount_apron_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "oceanic_plateau_edge_cell_fraction_ocean": (
            float(np.count_nonzero((marine.oceanic_plateau_edge_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "abyssal_hill_field_cell_fraction_ocean": (
            float(np.count_nonzero((marine.abyssal_hill_field_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "fracture_zone_rank_cell_fraction_ocean": (
            float(np.count_nonzero((marine.fracture_zone_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "abyssal_plain_fabric_rank_cell_fraction_ocean": (
            float(np.count_nonzero(
                (marine.abyssal_plain_fabric_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "coastal_plain_cell_fraction_land": (
            float(np.count_nonzero((coastal.coastal_plain_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "coastal_cliff_cell_fraction_land": (
            float(np.count_nonzero((coastal.coastal_cliff_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "shoreface_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.shoreface_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "barrier_lagoon_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.barrier_lagoon_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "estuary_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.estuary_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "delta_distributary_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.delta_distributary_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "estuary_funnel_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.estuary_funnel_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "barrier_spit_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.barrier_spit_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "delta_mouth_bar_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.delta_mouth_bar_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "estuary_tidal_channel_cell_fraction_ocean": (
            float(np.count_nonzero((coastal.estuary_tidal_channel_rank > 0.0) & ocean)
                  / max(int(np.count_nonzero(ocean)), 1))
        ),
        "coastal_depositional_plain_rank_cell_fraction_land": (
            float(np.count_nonzero((coastal.coastal_depositional_plain_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "strandplain_rank_cell_fraction_land": (
            float(np.count_nonzero((coastal.strandplain_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "tidal_flat_rank_cell_fraction_land": (
            float(np.count_nonzero((coastal.tidal_flat_rank > 0.0) & land)
                  / max(int(np.count_nonzero(land)), 1))
        ),
        "ocean_coast_distance_p95_passes": pct(
            marine.ocean_coast_distance[marine.ocean_coast_distance >= 0.0], 95.0),
        "land_coast_distance_p95_passes": pct(
            micro.land_coast_distance[micro.land_coast_distance >= 0.0], 95.0),
        "parent_rel_p05_m": pct(parent_rel, 5.0),
        "parent_rel_p50_m": pct(parent_rel, 50.0),
        "parent_rel_p95_m": pct(parent_rel, 95.0),
        "refined_rel_p05_m": pct(refined_rel, 5.0),
        "refined_rel_p50_m": pct(refined_rel, 50.0),
        "refined_rel_p95_m": pct(refined_rel, 95.0),
        "process_area_fractions": {
            key: frac(mask) for key, mask in process_masks.items()
        },
        "parent_p110a_summary": parent_metrics.get(
            "p110a_modern_planform", {}).get("summary", {}),
        "parent_p110a_warning_flags": parent_metrics.get(
            "p110a_modern_planform", {}).get("warning_flags", []),
    }


def _p107_compatible_metrics(
    refined_metrics: dict[str, Any],
    parent_metrics: dict[str, Any],
    arrays_out: Path,
    target_arrays: dict[str, np.ndarray],
    fields: dict[str, str],
) -> dict[str, Any]:
    manifest_fields = dict(fields)
    manifest_fields.setdefault("terrain.elevation_m", "field__terrain_elevation_m")
    manifest_fields.setdefault("tectonics.plate_id", "field__tectonics_plate_id")
    manifest_fields.setdefault("crust.type", "field__crust_type")
    manifest_fields.setdefault("crust.age_myr", "field__crust_age_myr")
    p110a = parent_metrics.get("p110a_modern_planform", {})
    p109 = parent_metrics.get("p109_hypsometry_comparison", {})
    return {
        "schema": "aevum.p107_terminal_world_audit.v1",
        "derived_from_schema": SCHEMA,
        "cells": int(refined_metrics["target_cells"]),
        "sea_level_m": float(refined_metrics["sea_level_m"]),
        "array_archive": {
            "path": str(arrays_out),
            "manifest": {
                "fields": manifest_fields,
                "arrays": sorted(str(key) for key in target_arrays),
            },
        },
        "p109_hypsometry_comparison": p109,
        "p110a_modern_planform": p110a,
        "selected_snapshot_refinement": refined_metrics,
    }


def _render_refinement_assets(
    grid: SphereGrid,
    parent_rel: np.ndarray,
    refined_rel: np.ndarray,
    detail_delta: np.ndarray,
    amplitude: np.ndarray,
    marine: _MarineMicrogeomorphology,
    micro: _Microgeomorphology,
    coastal: _CoastalMorphology,
    outdir: Path,
    *,
    width: int,
    height: int,
    render_groups: frozenset[str] | tuple[str, ...] | None = None,
) -> dict[str, Path]:
    groups = _normalize_render_groups(tuple(render_groups or ("all",)))
    if "all" not in groups:
        return _render_refinement_assets_selected(
            grid,
            parent_rel,
            refined_rel,
            detail_delta,
            amplitude,
            marine,
            micro,
            coastal,
            outdir,
            width=width,
            height=height,
            render_groups=groups,
        )
    outdir.mkdir(parents=True, exist_ok=True)
    parent_r = render.to_raster_continuous(
        grid, parent_rel, width=width, height=height, preserve_sign=True)
    refined_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid, detail_delta, width=width, height=height)
    amp_r = render.to_raster_continuous(
        grid, amplitude, width=width, height=height)

    assets: dict[str, Path] = {}
    assets["parent_vs_refined_elevation.png"] = _render_parent_refined_delta(
        parent_r, refined_r, delta_r, outdir / "parent_vs_refined_elevation.png")
    assets["refinement_delta_m.png"] = _render_delta(
        delta_r, outdir / "refinement_delta_m.png")
    assets["refinement_amplitude_m.png"] = _render_amplitude(
        amp_r, outdir / "refinement_amplitude_m.png")
    assets["refinement_zoom_sheet.png"] = _render_zoom_sheet(
        grid,
        parent_rel,
        refined_rel,
        detail_delta,
        outdir / "refinement_zoom_sheet.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_hydrology.png"] = _render_hydrology(
        grid,
        refined_rel,
        micro,
        outdir / "selected_snapshot_hydrology.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_hydrology_zoom_sheet.png"] = _render_hydrology_zoom_sheet(
        grid,
        refined_rel,
        micro,
        outdir / "selected_snapshot_hydrology_zoom_sheet.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_drainage_basins.png"] = _render_drainage_basins(
        grid,
        refined_rel,
        micro,
        outdir / "selected_snapshot_drainage_basins.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_fluvial_lacustrine_microshapes.png"] = (
        _render_fluvial_lacustrine_microshapes(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_fluvial_lacustrine_microshapes.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png"] = (
        _render_fluvial_lacustrine_microshapes_zoom_sheet(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_fluvial_microrelief_delta.png"] = (
        _render_fluvial_microrelief_delta(
            grid,
            micro,
            outdir / "selected_snapshot_fluvial_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_fluvial_microrelief_zoom_sheet.png"] = (
        _render_fluvial_microrelief_zoom_sheet(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_fluvial_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_lowland_alluvial_microrelief_delta.png"] = (
        _render_lowland_alluvial_microrelief_delta(
            grid,
            micro,
            outdir / "selected_snapshot_lowland_alluvial_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png"] = (
        _render_lowland_alluvial_microrelief_zoom_sheet(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_marine_microgeomorphology.png"] = _render_marine_microgeomorphology(
        grid,
        refined_rel,
        marine,
        outdir / "selected_snapshot_marine_microgeomorphology.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_marine_zoom_sheet.png"] = _render_marine_zoom_sheet(
        grid,
        refined_rel,
        marine,
        outdir / "selected_snapshot_marine_zoom_sheet.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_marine_object_classes.png"] = _render_marine_object_classes(
        grid,
        refined_rel,
        marine,
        outdir / "selected_snapshot_marine_object_classes.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_shelf_slope_microrelief_delta.png"] = (
        _render_shelf_slope_microrelief_delta(
            grid,
            marine,
            outdir / "selected_snapshot_shelf_slope_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_shelf_slope_microrelief_zoom_sheet.png"] = (
        _render_shelf_slope_microrelief_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_shelf_slope_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_deep_ocean_fabric_delta.png"] = (
        _render_deep_ocean_fabric_delta(
            grid,
            marine,
            outdir / "selected_snapshot_deep_ocean_fabric_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_deep_ocean_fabric_zoom_sheet.png"] = (
        _render_deep_ocean_fabric_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_deep_ocean_fabric_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_submarine_highlands.png"] = (
        _render_submarine_highlands(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_submarine_highlands.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_submarine_highlands_zoom_sheet.png"] = (
        _render_submarine_highlands_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_submarine_highlands_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_island_atoll_candidates.png"] = (
        _render_island_atoll_candidates(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_island_atoll_candidates.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_process_island_promotion.png"] = (
        _render_process_island_promotion(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_process_island_promotion.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_process_island_promotion_zoom_sheet.png"] = (
        _render_process_island_promotion_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_process_island_promotion_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_island_atoll_microshapes.png"] = (
        _render_island_atoll_microshapes(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_island_atoll_microshapes.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_island_atoll_microshapes_zoom_sheet.png"] = (
        _render_island_atoll_microshapes_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_island_atoll_microshapes_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_island_atoll_microrelief_delta.png"] = (
        _render_island_atoll_microrelief_delta(
            grid,
            marine,
            outdir / "selected_snapshot_island_atoll_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_island_atoll_microrelief_zoom_sheet.png"] = (
        _render_island_atoll_microrelief_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_island_atoll_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_reef_atoll_morphology.png"] = (
        _render_reef_atoll_morphology(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_reef_atoll_morphology.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_reef_atoll_zoom_sheet.png"] = (
        _render_reef_atoll_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_reef_atoll_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_morphology.png"] = _render_coastal_morphology(
        grid,
        refined_rel,
        coastal,
        outdir / "selected_snapshot_coastal_morphology.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_coastal_zoom_sheet.png"] = _render_coastal_zoom_sheet(
        grid,
        refined_rel,
        coastal,
        outdir / "selected_snapshot_coastal_zoom_sheet.png",
        width=width,
        height=height,
    )
    assets["selected_snapshot_coastal_process_linework.png"] = (
        _render_coastal_process_linework(
            grid,
            refined_rel,
            coastal,
            micro,
            outdir / "selected_snapshot_coastal_process_linework.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_process_linework_zoom_sheet.png"] = (
        _render_coastal_process_linework_zoom_sheet(
            grid,
            refined_rel,
            coastal,
            micro,
            outdir / "selected_snapshot_coastal_process_linework_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_process_microrelief_delta.png"] = (
        _render_coastal_process_microrelief_delta(
            grid,
            coastal,
            outdir / "selected_snapshot_coastal_process_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_process_microrelief_zoom_sheet.png"] = (
        _render_coastal_process_microrelief_zoom_sheet(
            grid,
            refined_rel,
            coastal,
            micro,
            outdir / "selected_snapshot_coastal_process_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_depositional_microrelief_delta.png"] = (
        _render_coastal_depositional_microrelief_delta(
            grid,
            coastal,
            outdir / "selected_snapshot_coastal_depositional_microrelief_delta.png",
            width=width,
            height=height,
        )
    )
    assets["selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png"] = (
        _render_coastal_depositional_microrelief_zoom_sheet(
            grid,
            refined_rel,
            coastal,
            outdir / "selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png",
            width=width,
            height=height,
        )
    )
    return assets


def _render_refinement_assets_selected(
    grid: SphereGrid,
    parent_rel: np.ndarray,
    refined_rel: np.ndarray,
    detail_delta: np.ndarray,
    amplitude: np.ndarray,
    marine: _MarineMicrogeomorphology,
    micro: _Microgeomorphology,
    coastal: _CoastalMorphology,
    outdir: Path,
    *,
    width: int,
    height: int,
    render_groups: frozenset[str],
) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    assets: dict[str, Path] = {}

    if _render_group_enabled(render_groups, "base"):
        parent_r = render.to_raster_continuous(
            grid, parent_rel, width=width, height=height, preserve_sign=True)
        refined_r = render.to_raster_continuous(
            grid, refined_rel, width=width, height=height, preserve_sign=True)
        delta_r = render.to_raster_continuous(
            grid, detail_delta, width=width, height=height)
        amp_r = render.to_raster_continuous(
            grid, amplitude, width=width, height=height)
        assets["parent_vs_refined_elevation.png"] = _render_parent_refined_delta(
            parent_r,
            refined_r,
            delta_r,
            outdir / "parent_vs_refined_elevation.png",
        )
        assets["refinement_delta_m.png"] = _render_delta(
            delta_r, outdir / "refinement_delta_m.png")
        assets["refinement_amplitude_m.png"] = _render_amplitude(
            amp_r, outdir / "refinement_amplitude_m.png")
        assets["refinement_zoom_sheet.png"] = _render_zoom_sheet(
            grid,
            parent_rel,
            refined_rel,
            detail_delta,
            outdir / "refinement_zoom_sheet.png",
            width=width,
            height=height,
        )

    if _render_group_enabled(render_groups, "hydrology"):
        assets["selected_snapshot_hydrology.png"] = _render_hydrology(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_hydrology.png",
            width=width,
            height=height,
        )
        assets["selected_snapshot_hydrology_zoom_sheet.png"] = (
            _render_hydrology_zoom_sheet(
                grid,
                refined_rel,
                micro,
                outdir / "selected_snapshot_hydrology_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_drainage_basins.png"] = _render_drainage_basins(
            grid,
            refined_rel,
            micro,
            outdir / "selected_snapshot_drainage_basins.png",
            width=width,
            height=height,
        )
        assets["selected_snapshot_fluvial_lacustrine_microshapes.png"] = (
            _render_fluvial_lacustrine_microshapes(
                grid,
                refined_rel,
                micro,
                outdir / "selected_snapshot_fluvial_lacustrine_microshapes.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png"] = (
            _render_fluvial_lacustrine_microshapes_zoom_sheet(
                grid,
                refined_rel,
                micro,
                outdir / "selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_fluvial_microrelief_delta.png"] = (
            _render_fluvial_microrelief_delta(
                grid,
                micro,
                outdir / "selected_snapshot_fluvial_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_fluvial_microrelief_zoom_sheet.png"] = (
            _render_fluvial_microrelief_zoom_sheet(
                grid,
                refined_rel,
                micro,
                outdir / "selected_snapshot_fluvial_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_lowland_alluvial_microrelief_delta.png"] = (
            _render_lowland_alluvial_microrelief_delta(
                grid,
                micro,
                outdir / "selected_snapshot_lowland_alluvial_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png"] = (
            _render_lowland_alluvial_microrelief_zoom_sheet(
                grid,
                refined_rel,
                micro,
                outdir / "selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "marine"):
        assets["selected_snapshot_marine_microgeomorphology.png"] = (
            _render_marine_microgeomorphology(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_marine_microgeomorphology.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_marine_zoom_sheet.png"] = _render_marine_zoom_sheet(
            grid,
            refined_rel,
            marine,
            outdir / "selected_snapshot_marine_zoom_sheet.png",
            width=width,
            height=height,
        )
        assets["selected_snapshot_marine_object_classes.png"] = (
            _render_marine_object_classes(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_marine_object_classes.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "shelf"):
        assets["selected_snapshot_shelf_slope_microrelief_delta.png"] = (
            _render_shelf_slope_microrelief_delta(
                grid,
                marine,
                outdir / "selected_snapshot_shelf_slope_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_shelf_slope_microrelief_zoom_sheet.png"] = (
            _render_shelf_slope_microrelief_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_shelf_slope_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "deep-ocean"):
        assets["selected_snapshot_deep_ocean_fabric_delta.png"] = (
            _render_deep_ocean_fabric_delta(
                grid,
                marine,
                outdir / "selected_snapshot_deep_ocean_fabric_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_deep_ocean_fabric_zoom_sheet.png"] = (
            _render_deep_ocean_fabric_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_deep_ocean_fabric_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "submarine"):
        assets["selected_snapshot_submarine_highlands.png"] = (
            _render_submarine_highlands(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_submarine_highlands.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_submarine_highlands_zoom_sheet.png"] = (
            _render_submarine_highlands_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_submarine_highlands_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "island-atoll"):
        assets["selected_snapshot_island_atoll_candidates.png"] = (
            _render_island_atoll_candidates(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_island_atoll_candidates.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_process_island_promotion.png"] = (
            _render_process_island_promotion(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_process_island_promotion.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_process_island_promotion_zoom_sheet.png"] = (
            _render_process_island_promotion_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_process_island_promotion_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_island_atoll_microshapes.png"] = (
            _render_island_atoll_microshapes(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_island_atoll_microshapes.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_island_atoll_microshapes_zoom_sheet.png"] = (
            _render_island_atoll_microshapes_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_island_atoll_microshapes_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_island_atoll_microrelief_delta.png"] = (
            _render_island_atoll_microrelief_delta(
                grid,
                marine,
                outdir / "selected_snapshot_island_atoll_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_island_atoll_microrelief_zoom_sheet.png"] = (
            _render_island_atoll_microrelief_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_island_atoll_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_reef_atoll_morphology.png"] = (
            _render_reef_atoll_morphology(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_reef_atoll_morphology.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_reef_atoll_zoom_sheet.png"] = (
            _render_reef_atoll_zoom_sheet(
                grid,
                refined_rel,
                marine,
                outdir / "selected_snapshot_reef_atoll_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    if _render_group_enabled(render_groups, "coastal"):
        assets["selected_snapshot_coastal_morphology.png"] = (
            _render_coastal_morphology(
                grid,
                refined_rel,
                coastal,
                outdir / "selected_snapshot_coastal_morphology.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_zoom_sheet.png"] = (
            _render_coastal_zoom_sheet(
                grid,
                refined_rel,
                coastal,
                outdir / "selected_snapshot_coastal_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_process_linework.png"] = (
            _render_coastal_process_linework(
                grid,
                refined_rel,
                coastal,
                micro,
                outdir / "selected_snapshot_coastal_process_linework.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_process_linework_zoom_sheet.png"] = (
            _render_coastal_process_linework_zoom_sheet(
                grid,
                refined_rel,
                coastal,
                micro,
                outdir / "selected_snapshot_coastal_process_linework_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_process_microrelief_delta.png"] = (
            _render_coastal_process_microrelief_delta(
                grid,
                coastal,
                outdir / "selected_snapshot_coastal_process_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_process_microrelief_zoom_sheet.png"] = (
            _render_coastal_process_microrelief_zoom_sheet(
                grid,
                refined_rel,
                coastal,
                micro,
                outdir / "selected_snapshot_coastal_process_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_depositional_microrelief_delta.png"] = (
            _render_coastal_depositional_microrelief_delta(
                grid,
                coastal,
                outdir / "selected_snapshot_coastal_depositional_microrelief_delta.png",
                width=width,
                height=height,
            )
        )
        assets["selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png"] = (
            _render_coastal_depositional_microrelief_zoom_sheet(
                grid,
                refined_rel,
                coastal,
                outdir / "selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png",
                width=width,
                height=height,
            )
        )

    return assets


def _render_parent_refined_delta(
    parent_r: np.ndarray,
    refined_r: np.ndarray,
    delta_r: np.ndarray,
    path: Path,
) -> Path:
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), constrained_layout=True)
    im0 = render.render_elevation_raster(axes[0], parent_r, title="Parent elevation")
    render.add_elevation_colorbar(fig, axes[0], im0)
    im1 = render.render_elevation_raster(axes[1], refined_r, title="Refined elevation")
    render.add_elevation_colorbar(fig, axes[1], im1)
    vmax = max(float(np.nanpercentile(np.abs(delta_r), 98.0)), 50.0)
    im2 = axes[2].imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    axes[2].set_title("72000 selected-snapshot detail delta (m)")
    fig.colorbar(im2, ax=axes[2], shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_delta(delta_r: np.ndarray, path: Path) -> Path:
    vmax = max(float(np.nanpercentile(np.abs(delta_r), 98.0)), 50.0)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 selected-snapshot detail delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_fluvial_microrelief_delta(
    grid: SphereGrid,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        micro.fluvial_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 fluvial / lacustrine microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_fluvial_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        micro.fluvial_microrelief_delta,
        width=width,
        height=height,
    )
    centers = _fluvial_lacustrine_zoom_centers(grid, micro)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=34.0, lat_height=20.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=34.0, lat_height=20.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        _draw_fluvial_lacustrine_microshapes(
            axes[row, 1],
            grid,
            micro,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: applied microrelief delta")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_lowland_alluvial_microrelief_delta(
    grid: SphereGrid,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        micro.lowland_alluvial_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        18.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 lowland / alluvial microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_lowland_alluvial_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        micro.lowland_alluvial_microrelief_delta,
        width=width,
        height=height,
    )
    fan_r = render.to_raster_continuous(
        grid, micro.alluvial_fan_rank, width=width, height=height)
    plain_r = render.to_raster_continuous(
        grid, micro.lowland_plain_rank, width=width, height=height)
    piedmont_r = render.to_raster_continuous(
        grid, micro.piedmont_apron_rank, width=width, height=height)
    centers = _lowland_alluvial_zoom_centers(grid, micro)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        18.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=38.0, lat_height=22.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=38.0, lat_height=22.0)
        fan_crop, _ = _crop_lon_lat(
            fan_r, lon, lat, lon_width=38.0, lat_height=22.0)
        plain_crop, _ = _crop_lon_lat(
            plain_r, lon, lat, lon_width=38.0, lat_height=22.0)
        piedmont_crop, _ = _crop_lon_lat(
            piedmont_r, lon, lat, lon_width=38.0, lat_height=22.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        for crop, cmap, threshold, alpha_scale in (
            (plain_crop, "Greens", 0.10, 0.34),
            (piedmont_crop, "Purples", 0.10, 0.42),
            (fan_crop, "YlOrBr", 0.10, 0.52),
        ):
            masked = np.ma.masked_where(crop <= threshold, crop)
            axes[row, 1].imshow(
                masked,
                cmap=cmap,
                vmin=threshold,
                vmax=1.0,
                alpha=np.clip(crop * alpha_scale, 0.0, min(0.88, alpha_scale + 0.12)),
                extent=extent,
            )
        axes[row, 1].set_title(f"{label}: alluvial delta + ranks")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_amplitude(amp_r: np.ndarray, path: Path) -> Path:
    vmax = max(float(np.nanpercentile(amp_r, 98.0)), 50.0)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        amp_r,
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 selected-snapshot process detail amplitude (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_drainage_basins(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    basin_r = render.to_raster(
        grid,
        micro.drainage_basin_id.astype(np.float64),
        width=width,
        height=height,
    )
    floodplain_r = render.to_raster_continuous(
        grid, micro.floodplain_rank, width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax, elev_r, title="72000 selected-snapshot drainage basins")
    basin = np.ma.masked_where(basin_r <= 0.0, basin_r)
    ax.imshow(
        basin,
        cmap="tab20",
        alpha=np.where(basin_r > 0.0, 0.28, 0.0),
        extent=[-180, 180, -90, 90],
        zorder=3,
    )
    floodplain = np.ma.masked_where(floodplain_r <= 0.05, floodplain_r)
    ax.imshow(
        floodplain,
        cmap="Greens",
        vmin=0.05,
        vmax=1.0,
        alpha=np.clip(floodplain_r * 0.46, 0.0, 0.48),
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    _draw_river_paths(
        ax,
        grid,
        micro,
        min_rank=0.20,
        rank_values=micro.basin_trunk_rank,
        meander_values=micro.meander_belt_rank,
        color="#063d8f",
        linewidth_base=0.35,
        linewidth_scale=1.70,
        alpha_base=0.48,
        alpha_scale=0.44,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_marine_microgeomorphology(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    shelf_r = render.to_raster_continuous(
        grid, marine.shelf_break_rank, width=width, height=height)
    reef_r = render.to_raster(
        grid, marine.reef_atoll_mask.astype(np.float64), width=width, height=height)
    shoal_r = render.to_raster(
        grid, marine.marine_shoal_mask.astype(np.float64), width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax, elev_r, title="72000 selected-snapshot marine microgeomorphology")
    shelf = np.ma.masked_where(shelf_r < 0.68, shelf_r)
    ax.imshow(
        shelf,
        cmap="Greys",
        vmin=0.68,
        vmax=1.0,
        alpha=np.clip((shelf_r - 0.55) * 0.58, 0.0, 0.38),
        extent=[-180, 180, -90, 90],
        zorder=3,
    )
    shoal = np.ma.masked_where(shoal_r <= 0.5, shoal_r)
    ax.imshow(
        shoal,
        cmap="Wistia",
        vmin=0.0,
        vmax=1.0,
        alpha=0.48,
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    reef = np.ma.masked_where(reef_r <= 0.5, reef_r)
    ax.imshow(
        reef,
        cmap="winter",
        vmin=0.0,
        vmax=1.0,
        alpha=0.92,
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_marine_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    shelf_r = render.to_raster_continuous(
        grid, marine.shelf_break_rank, width=width, height=height)
    reef_r = render.to_raster(
        grid, marine.reef_atoll_mask.astype(np.float64), width=width, height=height)
    shoal_r = render.to_raster(
        grid, marine.marine_shoal_mask.astype(np.float64), width=width, height=height)
    island_candidate_r = render.to_raster_continuous(
        grid, marine.island_candidate_rank, width=width, height=height)
    atoll_candidate_r = render.to_raster_continuous(
        grid, marine.atoll_candidate_rank, width=width, height=height)
    reef_rim_r = render.to_raster_continuous(
        grid, marine.reef_rim_rank, width=width, height=height)
    lagoon_r = render.to_raster_continuous(
        grid, marine.atoll_lagoon_rank, width=width, height=height)
    fringing_r = render.to_raster_continuous(
        grid, marine.fringing_reef_rank, width=width, height=height)
    centers = _marine_zoom_centers(grid, refined_rel, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.9 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(elev_r, lon, lat, lon_width=54.0, lat_height=28.0)
        shelf_crop, _ = _crop_lon_lat(shelf_r, lon, lat, lon_width=54.0, lat_height=28.0)
        reef_crop, _ = _crop_lon_lat(reef_r, lon, lat, lon_width=54.0, lat_height=28.0)
        shoal_crop, _ = _crop_lon_lat(shoal_r, lon, lat, lon_width=54.0, lat_height=28.0)
        island_crop, _ = _crop_lon_lat(
            island_candidate_r, lon, lat, lon_width=54.0, lat_height=28.0)
        atoll_crop, _ = _crop_lon_lat(
            atoll_candidate_r, lon, lat, lon_width=54.0, lat_height=28.0)
        rim_crop, _ = _crop_lon_lat(
            reef_rim_r, lon, lat, lon_width=54.0, lat_height=28.0)
        lagoon_crop, _ = _crop_lon_lat(
            lagoon_r, lon, lat, lon_width=54.0, lat_height=28.0)
        fringing_crop, _ = _crop_lon_lat(
            fringing_r, lon, lat, lon_width=54.0, lat_height=28.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        shelf = np.ma.masked_where(shelf_crop < 0.68, shelf_crop)
        axes[row, 1].imshow(
            shelf,
            cmap="Greys",
            vmin=0.68,
            vmax=1.0,
            alpha=np.clip((shelf_crop - 0.55) * 0.62, 0.0, 0.42),
            extent=extent,
        )
        shoal = np.ma.masked_where(shoal_crop <= 0.5, shoal_crop)
        axes[row, 1].imshow(
            shoal,
            cmap="Wistia",
            vmin=0.0,
            vmax=1.0,
            alpha=0.50,
            extent=extent,
        )
        reef = np.ma.masked_where(reef_crop <= 0.5, reef_crop)
        axes[row, 1].imshow(
            reef,
            cmap="winter",
            vmin=0.0,
            vmax=1.0,
            alpha=0.92,
            extent=extent,
        )
        island = np.ma.masked_where(island_crop <= 0.12, island_crop)
        axes[row, 1].imshow(
            island,
            cmap="autumn",
            vmin=0.12,
            vmax=1.0,
            alpha=np.clip(island_crop * 0.66, 0.0, 0.72),
            extent=extent,
        )
        atoll = np.ma.masked_where(atoll_crop <= 0.12, atoll_crop)
        axes[row, 1].imshow(
            atoll,
            cmap="winter",
            vmin=0.12,
            vmax=1.0,
            alpha=np.clip(atoll_crop * 0.72, 0.0, 0.82),
            extent=extent,
        )
        lagoon = np.ma.masked_where(lagoon_crop <= 0.08, lagoon_crop)
        axes[row, 1].imshow(
            lagoon,
            cmap="Blues",
            vmin=0.08,
            vmax=1.0,
            alpha=np.clip(lagoon_crop * 0.56, 0.0, 0.64),
            extent=extent,
        )
        fringing = np.ma.masked_where(fringing_crop <= 0.08, fringing_crop)
        axes[row, 1].imshow(
            fringing,
            cmap="summer",
            vmin=0.08,
            vmax=1.0,
            alpha=np.clip(fringing_crop * 0.56, 0.0, 0.68),
            extent=extent,
        )
        rim = np.ma.masked_where(rim_crop <= 0.08, rim_crop)
        axes[row, 1].imshow(
            rim,
            cmap="Wistia",
            vmin=0.08,
            vmax=1.0,
            alpha=np.clip(rim_crop * 0.66, 0.0, 0.78),
            extent=extent,
        )
        axes[row, 1].set_title(f"{label}: shelf/shoal/reef/candidates")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.75)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_marine_object_classes(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    class_field = _marine_object_class_field(marine)
    class_r = render.to_raster(grid, class_field, width=width, height=height)
    masked = np.ma.masked_where(class_r <= 0.0, class_r - 1.0)
    cmap = matplotlib.colors.ListedColormap([
        "#5b6470",  # shelf break
        "#d18a22",  # seamount shoal
        "#b7a34a",  # oceanic plateau shoal
        "#7bbf6a",  # microcontinent shoal
        "#c95b3f",  # island arc shoal
        "#28d7b3",  # reef / atoll
    ])
    norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5, 6.5, 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax, elev_r, title="72000 selected-snapshot marine object classes")
    im = ax.imshow(
        masked,
        cmap=cmap,
        norm=norm,
        alpha=0.82,
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.72, ticks=np.arange(0, 6))
    cb.ax.set_yticklabels([
        "shelf break",
        "seamount",
        "plateau",
        "microcontinent",
        "island arc",
        "reef/atoll",
    ])
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _marine_object_class_field(marine: _MarineMicrogeomorphology) -> np.ndarray:
    out = np.zeros(marine.shelf_break_rank.shape, dtype=np.int16)
    out[marine.shelf_break_rank >= 0.68] = 1
    out[marine.seamount_shoal_mask] = 2
    out[marine.oceanic_plateau_shoal_mask] = 3
    out[marine.microcontinent_shoal_mask] = 4
    out[marine.island_arc_shoal_mask] = 5
    out[marine.reef_atoll_mask] = 6
    return out


def _render_shelf_slope_microrelief_delta(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        marine.shelf_slope_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 shelf / upper-slope microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_shelf_slope_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        marine.shelf_slope_microrelief_delta,
        width=width,
        height=height,
    )
    shelf_r = render.to_raster_continuous(
        grid, marine.shelf_break_rank, width=width, height=height)
    centers = _shelf_slope_microrelief_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=42.0, lat_height=22.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=42.0, lat_height=22.0)
        shelf_crop, _ = _crop_lon_lat(
            shelf_r, lon, lat, lon_width=42.0, lat_height=22.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        shelf = np.ma.masked_where(shelf_crop < 0.62, shelf_crop)
        axes[row, 1].imshow(
            shelf,
            cmap="Greys",
            vmin=0.62,
            vmax=1.0,
            alpha=np.clip((shelf_crop - 0.52) * 0.56, 0.0, 0.40),
            extent=extent,
        )
        axes[row, 1].set_title(f"{label}: applied shelf/slope delta")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_deep_ocean_fabric_delta(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        marine.deep_ocean_fabric_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        16.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 deep-ocean fabric delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_deep_ocean_fabric_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        marine.deep_ocean_fabric_delta,
        width=width,
        height=height,
    )
    fracture_r = render.to_raster_continuous(
        grid, marine.fracture_zone_rank, width=width, height=height)
    plain_r = render.to_raster_continuous(
        grid, marine.abyssal_plain_fabric_rank, width=width, height=height)
    centers = _deep_ocean_fabric_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        16.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=48.0, lat_height=26.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=48.0, lat_height=26.0)
        fracture_crop, _ = _crop_lon_lat(
            fracture_r, lon, lat, lon_width=48.0, lat_height=26.0)
        plain_crop, _ = _crop_lon_lat(
            plain_r, lon, lat, lon_width=48.0, lat_height=26.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: bathymetry")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        plain = np.ma.masked_where(plain_crop <= 0.10, plain_crop)
        axes[row, 1].imshow(
            plain,
            cmap="Greys",
            vmin=0.10,
            vmax=1.0,
            alpha=np.clip(plain_crop * 0.28, 0.0, 0.34),
            extent=extent,
            zorder=4,
        )
        fracture = np.ma.masked_where(fracture_crop <= 0.10, fracture_crop)
        axes[row, 1].imshow(
            fracture,
            cmap="magma",
            vmin=0.10,
            vmax=1.0,
            alpha=np.clip(fracture_crop * 0.44, 0.0, 0.52),
            extent=extent,
            zorder=5,
        )
        _draw_deep_ocean_fracture_linework(
            axes[row, 1],
            grid,
            marine,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: fabric delta + ranks")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _draw_deep_ocean_fracture_linework(
    ax,
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    *,
    extent: list[float] | None,
    zoom: bool,
) -> None:
    rank = np.asarray(marine.fracture_zone_rank, dtype=np.float64)
    active = rank > 0.12
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        active &= (
            (grid.lon >= lon_min)
            & (grid.lon <= lon_max)
            & (grid.lat >= lat_min)
            & (grid.lat <= lat_max)
        )
    if not np.any(active):
        return
    seen = np.zeros(grid.n, dtype=bool)
    components: list[tuple[float, np.ndarray]] = []
    for start in np.flatnonzero(active):
        if seen[int(start)]:
            continue
        stack = [int(start)]
        seen[int(start)] = True
        cells: list[int] = []
        while stack:
            cell = stack.pop()
            cells.append(cell)
            for nb in grid.neighbors[cell]:
                nb = int(nb)
                if active[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        if len(cells) < 3:
            continue
        comp = np.asarray(cells, dtype=np.int64)
        score = float(np.nanmax(rank[comp]) * np.sqrt(comp.size))
        components.append((score, comp))
    if not components:
        return
    components.sort(key=lambda item: item[0], reverse=True)
    max_lines = 10 if zoom else 42
    for _, comp in components[:max_lines]:
        lon = np.asarray(grid.lon[comp], dtype=np.float64)
        lat = np.asarray(grid.lat[comp], dtype=np.float64)
        center_lon = float(lon[np.nanargmax(rank[comp])])
        lon = center_lon + ((lon - center_lon + 180.0) % 360.0 - 180.0)
        if np.nanmax(lon) - np.nanmin(lon) > 36.0:
            continue
        pts = np.column_stack([lon, lat])
        mean = np.nanmean(pts, axis=0)
        centered = pts - mean
        if not np.all(np.isfinite(centered)):
            continue
        try:
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        axis = vh[0]
        order = np.argsort(centered @ axis)
        x = lon[order]
        y = lat[order]
        if x.size >= 5:
            kernel = np.ones(3, dtype=np.float64) / 3.0
            x = np.convolve(x, kernel, mode="same")
            y = np.convolve(y, kernel, mode="same")
            x[0] = lon[order][0]
            y[0] = lat[order][0]
            x[-1] = lon[order][-1]
            y[-1] = lat[order][-1]
        strength = float(np.clip(np.nanmax(rank[comp]), 0.0, 1.0))
        ax.plot(
            x,
            y,
            color="#25306f",
            linewidth=(0.56 if zoom else 0.34) + (1.05 if zoom else 0.72) * strength,
            alpha=0.46 + 0.32 * strength,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=9,
        )


def _render_submarine_highlands(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    peak_r = render.to_raster_continuous(
        grid, marine.seamount_peak_rank, width=width, height=height)
    apron_r = render.to_raster_continuous(
        grid, marine.seamount_apron_rank, width=width, height=height)
    plateau_edge_r = render.to_raster_continuous(
        grid, marine.oceanic_plateau_edge_rank, width=width, height=height)
    abyssal_r = render.to_raster_continuous(
        grid, marine.abyssal_hill_field_rank, width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 submarine highland ranks (topology preserved)",
    )
    abyssal = np.ma.masked_where(abyssal_r <= 0.08, abyssal_r)
    ax.imshow(
        abyssal,
        cmap="Greys",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(abyssal_r * 0.34, 0.0, 0.38),
        extent=[-180, 180, -90, 90],
        zorder=3,
    )
    apron = np.ma.masked_where(apron_r <= 0.08, apron_r)
    ax.imshow(
        apron,
        cmap="YlOrBr",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(apron_r * 0.44, 0.0, 0.50),
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    plateau = np.ma.masked_where(plateau_edge_r <= 0.08, plateau_edge_r)
    ax.imshow(
        plateau,
        cmap="Purples",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(plateau_edge_r * 0.56, 0.0, 0.64),
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    peak = np.ma.masked_where(peak_r <= 0.08, peak_r)
    ax.imshow(
        peak,
        cmap="autumn",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(peak_r * 0.72, 0.0, 0.82),
        extent=[-180, 180, -90, 90],
        zorder=6,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_submarine_highlands_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    peak_r = render.to_raster_continuous(
        grid, marine.seamount_peak_rank, width=width, height=height)
    apron_r = render.to_raster_continuous(
        grid, marine.seamount_apron_rank, width=width, height=height)
    plateau_edge_r = render.to_raster_continuous(
        grid, marine.oceanic_plateau_edge_rank, width=width, height=height)
    abyssal_r = render.to_raster_continuous(
        grid, marine.abyssal_hill_field_rank, width=width, height=height)
    centers = _submarine_highlands_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.55 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=44.0, lat_height=24.0)
        peak_crop, _ = _crop_lon_lat(
            peak_r, lon, lat, lon_width=44.0, lat_height=24.0)
        apron_crop, _ = _crop_lon_lat(
            apron_r, lon, lat, lon_width=44.0, lat_height=24.0)
        plateau_crop, _ = _crop_lon_lat(
            plateau_edge_r, lon, lat, lon_width=44.0, lat_height=24.0)
        abyssal_crop, _ = _crop_lon_lat(
            abyssal_r, lon, lat, lon_width=44.0, lat_height=24.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: bathymetry")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        for crop, cmap, threshold, alpha_scale in (
            (abyssal_crop, "Greys", 0.08, 0.50),
            (apron_crop, "YlOrBr", 0.08, 0.62),
            (plateau_crop, "Purples", 0.08, 0.70),
            (peak_crop, "autumn", 0.08, 0.86),
        ):
            masked = np.ma.masked_where(crop <= threshold, crop)
            axes[row, 1].imshow(
                masked,
                cmap=cmap,
                vmin=threshold,
                vmax=1.0,
                alpha=np.clip(crop * alpha_scale, 0.0, min(0.94, alpha_scale + 0.10)),
                extent=extent,
            )
        axes[row, 1].set_title(f"{label}: peak/apron/plateau/abyssal hill")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_island_atoll_candidates(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    island_r = render.to_raster_continuous(
        grid, marine.island_candidate_rank, width=width, height=height)
    atoll_r = render.to_raster_continuous(
        grid, marine.atoll_candidate_rank, width=width, height=height)
    island_core_r = render.to_raster(
        grid,
        marine.process_island_candidate_mask.astype(np.float64),
        width=width,
        height=height,
    )
    atoll_core_r = render.to_raster(
        grid,
        marine.atoll_candidate_mask.astype(np.float64),
        width=width,
        height=height,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 process island / atoll candidates (topology preserved)",
    )
    island = np.ma.masked_where(island_r <= 0.12, island_r)
    ax.imshow(
        island,
        cmap="autumn",
        vmin=0.12,
        vmax=1.0,
        alpha=np.clip(island_r * 0.72, 0.0, 0.76),
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    atoll = np.ma.masked_where(atoll_r <= 0.12, atoll_r)
    ax.imshow(
        atoll,
        cmap="winter",
        vmin=0.12,
        vmax=1.0,
        alpha=np.clip(atoll_r * 0.76, 0.0, 0.84),
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    core = np.ma.masked_where(island_core_r <= 0.5, island_core_r)
    ax.imshow(
        core,
        cmap=matplotlib.colors.ListedColormap(["#fff2a8"]),
        vmin=0.0,
        vmax=1.0,
        alpha=0.94,
        extent=[-180, 180, -90, 90],
        zorder=6,
    )
    atoll_core = np.ma.masked_where(atoll_core_r <= 0.5, atoll_core_r)
    ax.imshow(
        atoll_core,
        cmap=matplotlib.colors.ListedColormap(["#39ffd2"]),
        vmin=0.0,
        vmax=1.0,
        alpha=0.90,
        extent=[-180, 180, -90, 90],
        zorder=7,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_process_island_promotion(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    promotion_r = render.to_raster_continuous(
        grid,
        marine.process_island_promotion_rank,
        width=width,
        height=height,
    )
    promoted_r = render.to_raster(
        grid,
        marine.process_island_promoted_mask.astype(np.float64),
        width=width,
        height=height,
    )
    atoll_r = render.to_raster(
        grid,
        marine.atoll_islet_promoted_mask.astype(np.float64),
        width=width,
        height=height,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 opt-in process island promotion",
    )
    promotion = np.ma.masked_where(promotion_r <= 0.08, promotion_r)
    ax.imshow(
        promotion,
        cmap="autumn",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(promotion_r * 0.66, 0.0, 0.76),
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    promoted = np.ma.masked_where(promoted_r <= 0.5, promoted_r)
    ax.imshow(
        promoted,
        cmap=matplotlib.colors.ListedColormap(["#fff1a6"]),
        vmin=0.0,
        vmax=1.0,
        alpha=0.96,
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    atoll = np.ma.masked_where(atoll_r <= 0.5, atoll_r)
    ax.imshow(
        atoll,
        cmap=matplotlib.colors.ListedColormap(["#41f7d0"]),
        vmin=0.0,
        vmax=1.0,
        alpha=0.94,
        extent=[-180, 180, -90, 90],
        zorder=6,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_process_island_promotion_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    candidate_r = render.to_raster_continuous(
        grid, marine.island_candidate_rank, width=width, height=height)
    atoll_candidate_r = render.to_raster_continuous(
        grid, marine.atoll_candidate_rank, width=width, height=height)
    promotion_r = render.to_raster_continuous(
        grid,
        marine.process_island_promotion_rank,
        width=width,
        height=height,
    )
    promoted_r = render.to_raster(
        grid,
        marine.process_island_promoted_mask.astype(np.float64),
        width=width,
        height=height,
    )
    atoll_r = render.to_raster(
        grid,
        marine.atoll_islet_promoted_mask.astype(np.float64),
        width=width,
        height=height,
    )
    centers = _process_island_promotion_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.55 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=32.0, lat_height=18.0)
        candidate_crop, _ = _crop_lon_lat(
            candidate_r, lon, lat, lon_width=32.0, lat_height=18.0)
        atoll_candidate_crop, _ = _crop_lon_lat(
            atoll_candidate_r, lon, lat, lon_width=32.0, lat_height=18.0)
        promotion_crop, _ = _crop_lon_lat(
            promotion_r, lon, lat, lon_width=32.0, lat_height=18.0)
        promoted_crop, _ = _crop_lon_lat(
            promoted_r, lon, lat, lon_width=32.0, lat_height=18.0)
        atoll_crop, _ = _crop_lon_lat(
            atoll_r, lon, lat, lon_width=32.0, lat_height=18.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        for crop, cmap, threshold, alpha_scale in (
            (candidate_crop, "autumn", 0.12, 0.52),
            (atoll_candidate_crop, "winter", 0.12, 0.56),
            (promotion_crop, "Oranges", 0.08, 0.70),
        ):
            masked = np.ma.masked_where(crop <= threshold, crop)
            axes[row, 1].imshow(
                masked,
                cmap=cmap,
                vmin=threshold,
                vmax=1.0,
                alpha=np.clip(crop * alpha_scale, 0.0, min(0.92, alpha_scale + 0.10)),
                extent=extent,
            )
        promoted = np.ma.masked_where(promoted_crop <= 0.5, promoted_crop)
        axes[row, 1].imshow(
            promoted,
            cmap=matplotlib.colors.ListedColormap(["#fff1a6"]),
            vmin=0.0,
            vmax=1.0,
            alpha=0.98,
            extent=extent,
        )
        atoll = np.ma.masked_where(atoll_crop <= 0.5, atoll_crop)
        axes[row, 1].imshow(
            atoll,
            cmap=matplotlib.colors.ListedColormap(["#41f7d0"]),
            vmin=0.0,
            vmax=1.0,
            alpha=0.96,
            extent=extent,
        )
        axes[row, 1].set_title(f"{label}: candidates / promoted land")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_island_atoll_microshapes(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 island / atoll subcell microshapes",
    )
    _draw_island_atoll_microshapes(ax, grid, marine, zoom=False)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_island_atoll_microshapes_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    centers = _island_atoll_microshape_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.55 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=26.0, lat_height=15.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: raster elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        _draw_island_atoll_microshapes(
            axes[row, 1],
            grid,
            marine,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: subcell island / atoll symbols")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_island_atoll_microrelief_delta(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        marine.island_atoll_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 island / atoll microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_island_atoll_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        marine.island_atoll_microrelief_delta,
        width=width,
        height=height,
    )
    centers = _island_atoll_microshape_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.55 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=26.0, lat_height=15.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=26.0, lat_height=15.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        _draw_island_atoll_microshapes(
            axes[row, 1],
            grid,
            marine,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: applied island / atoll delta")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _draw_island_atoll_microshapes(
    ax,
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    *,
    extent: list[float] | None = None,
    zoom: bool,
) -> None:
    atoll_cells = _microshape_cells(
        grid,
        marine.atoll_microshape_rank,
        min_score=0.26,
        max_count=180 if zoom else 70,
        extent=extent,
        force=marine.atoll_islet_promoted_mask,
    )
    islet_cells = _microshape_cells(
        grid,
        marine.islet_microshape_rank,
        min_score=0.50,
        max_count=160 if zoom else 70,
        extent=extent,
        force=marine.process_island_promoted_mask & ~marine.atoll_islet_promoted_mask,
        exclude=marine.atoll_islet_promoted_mask,
    )
    for cell in atoll_cells:
        _draw_atoll_microshape(ax, grid, marine, int(cell), zoom=zoom)
    for cell in islet_cells:
        _draw_islet_microshape(ax, grid, marine, int(cell), zoom=zoom)


def _microshape_cells(
    grid: SphereGrid,
    rank: np.ndarray,
    *,
    min_score: float,
    max_count: int,
    extent: list[float] | None,
    force: np.ndarray,
    exclude: np.ndarray | None = None,
) -> np.ndarray:
    values = np.asarray(rank, dtype=np.float64)
    forced = np.asarray(force, dtype=bool)
    mask = (values >= float(min_score)) | forced
    if exclude is not None:
        mask &= ~np.asarray(exclude, dtype=bool)
    cells = np.flatnonzero(mask)
    if cells.size == 0:
        return cells
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        inside = (
            (grid.lon[cells] >= lon_min)
            & (grid.lon[cells] <= lon_max)
            & (grid.lat[cells] >= lat_min)
            & (grid.lat[cells] <= lat_max)
        )
        cells = cells[inside]
    if cells.size <= int(max_count):
        return cells
    forced_cells = cells[forced[cells]]
    remaining = cells[~forced[cells]]
    keep_count = max(0, int(max_count) - int(forced_cells.size))
    if remaining.size and keep_count:
        ranked = remaining[np.argsort(values[remaining])[-keep_count:]]
        return np.concatenate([forced_cells, ranked])
    return forced_cells[:int(max_count)]


def _draw_islet_microshape(
    ax,
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(marine.islet_microshape_rank[cell], 0.0, 1.0))
    promoted = bool(marine.process_island_promoted_mask[cell])
    angle = _microshape_angle(cell)
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    base = 0.42 if zoom else 0.58
    width = base * lon_scale * (0.90 + 0.58 * rank)
    height = base * (0.38 + 0.28 * rank)
    face = "#fff1a6" if promoted else "#ffe7a3"
    edge = "#1f5f42" if promoted else "#8c7a3c"
    alpha = 0.92 if promoted else 0.54
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        width,
        height,
        angle=angle,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.05 if zoom else 0.62,
        alpha=alpha,
        zorder=8,
    ))
    cap_width = width * 0.54
    cap_height = height * 0.38
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        cap_width,
        cap_height,
        angle=angle + 8.0,
        facecolor="#2f8c5e",
        edgecolor="none",
        alpha=0.42 if promoted else 0.22,
        zorder=9,
    ))


def _draw_atoll_microshape(
    ax,
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(marine.atoll_microshape_rank[cell], 0.0, 1.0))
    promoted = bool(marine.atoll_islet_promoted_mask[cell])
    angle = _microshape_angle(cell) + 18.0
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    base = 0.48 if zoom else 0.66
    outer_w = base * lon_scale * (1.08 + 0.56 * rank)
    outer_h = base * (0.70 + 0.30 * rank)
    lagoon_w = outer_w * 0.54
    lagoon_h = outer_h * 0.46
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        outer_w,
        outer_h,
        angle=angle,
        facecolor=(0.18, 0.94, 0.80, 0.14 if promoted else 0.09),
        edgecolor="#5df2d0",
        linewidth=1.20 if zoom else 0.72,
        alpha=0.92 if promoted else 0.58,
        zorder=8,
    ))
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        lagoon_w,
        lagoon_h,
        angle=angle,
        facecolor="#66bde9",
        edgecolor="#2d8fb9",
        linewidth=0.55 if zoom else 0.32,
        alpha=0.48 if promoted else 0.30,
        zorder=9,
    ))
    theta = np.radians(angle)
    dx = 0.30 * outer_w * np.cos(theta)
    dy = 0.30 * outer_h * np.sin(theta)
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon + dx, lat + dy),
        outer_w * 0.34,
        outer_h * 0.22,
        angle=angle + 4.0,
        facecolor="#fff1a6",
        edgecolor="#7d6d3c",
        linewidth=0.70 if zoom else 0.42,
        alpha=0.95 if promoted else 0.62,
        zorder=10,
    ))


def _microshape_angle(cell: int) -> float:
    return float((int(cell) * 137.50776405003785) % 180.0)


def _render_reef_atoll_morphology(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    rim_r = render.to_raster_continuous(
        grid, marine.reef_rim_rank, width=width, height=height)
    lagoon_r = render.to_raster_continuous(
        grid, marine.atoll_lagoon_rank, width=width, height=height)
    fringing_r = render.to_raster_continuous(
        grid, marine.fringing_reef_rank, width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 reef / atoll morphology ranks (topology preserved)",
    )
    _plot_single_color_rank_overlay(
        ax,
        lagoon_r,
        threshold=0.12,
        color="#0b74a5",
        alpha_scale=0.72,
        alpha_max=0.78,
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    _draw_rank_point_overlay(
        ax,
        grid,
        marine.atoll_lagoon_rank,
        threshold=0.20,
        color="#0b74a5",
        edgecolor="#e7fbff",
        marker="o",
        extent=[-180, 180, -90, 90],
        max_count=220,
        size_base=7.0,
        size_scale=20.0,
        zorder=8,
    )
    _plot_single_color_rank_overlay(
        ax,
        fringing_r,
        threshold=0.12,
        color="#19a974",
        alpha_scale=0.70,
        alpha_max=0.78,
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    _draw_rank_point_overlay(
        ax,
        grid,
        marine.fringing_reef_rank,
        threshold=0.18,
        color="#19a974",
        edgecolor="#043d2a",
        marker="^",
        extent=[-180, 180, -90, 90],
        max_count=300,
        size_base=7.0,
        size_scale=18.0,
        zorder=9,
    )
    _plot_single_color_rank_overlay(
        ax,
        rim_r,
        threshold=0.14,
        color="#ffd23f",
        alpha_scale=0.86,
        alpha_max=0.94,
        extent=[-180, 180, -90, 90],
        zorder=6,
    )
    _draw_rank_point_overlay(
        ax,
        grid,
        marine.reef_rim_rank,
        threshold=0.20,
        color="#ffd23f",
        edgecolor="#6b4300",
        marker="o",
        extent=[-180, 180, -90, 90],
        max_count=320,
        size_base=8.0,
        size_scale=22.0,
        zorder=10,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_reef_atoll_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    rim_r = render.to_raster_continuous(
        grid, marine.reef_rim_rank, width=width, height=height)
    lagoon_r = render.to_raster_continuous(
        grid, marine.atoll_lagoon_rank, width=width, height=height)
    fringing_r = render.to_raster_continuous(
        grid, marine.fringing_reef_rank, width=width, height=height)
    centers = _reef_atoll_zoom_centers(grid, marine)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=32.0, lat_height=18.0)
        rim_crop, _ = _crop_lon_lat(
            rim_r, lon, lat, lon_width=32.0, lat_height=18.0)
        lagoon_crop, _ = _crop_lon_lat(
            lagoon_r, lon, lat, lon_width=32.0, lat_height=18.0)
        fringing_crop, _ = _crop_lon_lat(
            fringing_r, lon, lat, lon_width=32.0, lat_height=18.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        _plot_single_color_rank_overlay(
            axes[row, 1],
            lagoon_crop,
            threshold=0.12,
            color="#0b74a5",
            alpha_scale=0.86,
            alpha_max=0.92,
            extent=extent,
            zorder=4,
        )
        _draw_rank_point_overlay(
            axes[row, 1],
            grid,
            marine.atoll_lagoon_rank,
            threshold=0.20,
            color="#0b74a5",
            edgecolor="#e7fbff",
            marker="o",
            extent=extent,
            max_count=80,
            size_base=22.0,
            size_scale=42.0,
            zorder=8,
        )
        _plot_single_color_rank_overlay(
            axes[row, 1],
            fringing_crop,
            threshold=0.12,
            color="#19a974",
            alpha_scale=0.82,
            alpha_max=0.90,
            extent=extent,
            zorder=5,
        )
        _draw_rank_point_overlay(
            axes[row, 1],
            grid,
            marine.fringing_reef_rank,
            threshold=0.18,
            color="#19a974",
            edgecolor="#043d2a",
            marker="^",
            extent=extent,
            max_count=100,
            size_base=22.0,
            size_scale=38.0,
            zorder=9,
        )
        _plot_single_color_rank_overlay(
            axes[row, 1],
            rim_crop,
            threshold=0.14,
            color="#ffd23f",
            alpha_scale=0.98,
            alpha_max=0.98,
            extent=extent,
            zorder=6,
        )
        _draw_rank_point_overlay(
            axes[row, 1],
            grid,
            marine.reef_rim_rank,
            threshold=0.20,
            color="#ffd23f",
            edgecolor="#6b4300",
            marker="o",
            extent=extent,
            max_count=120,
            size_base=24.0,
            size_scale=44.0,
            zorder=10,
        )
        axes[row, 1].set_title(f"{label}: lagoon/fringing/rim")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _plot_single_color_rank_overlay(
    ax,
    rank_r: np.ndarray,
    *,
    threshold: float,
    color: str,
    alpha_scale: float,
    alpha_max: float,
    extent: list[float],
    zorder: int,
) -> None:
    masked = np.ma.masked_where(rank_r <= float(threshold), rank_r)
    ax.imshow(
        masked,
        cmap=matplotlib.colors.ListedColormap([color]),
        vmin=float(threshold),
        vmax=1.0,
        alpha=np.clip(rank_r * float(alpha_scale), 0.0, float(alpha_max)),
        extent=extent,
        zorder=zorder,
    )


def _draw_rank_point_overlay(
    ax,
    grid: SphereGrid,
    rank: np.ndarray,
    *,
    threshold: float,
    color: str,
    edgecolor: str,
    marker: str,
    extent: list[float],
    max_count: int,
    size_base: float,
    size_scale: float,
    zorder: int,
) -> None:
    values = np.asarray(rank, dtype=np.float64)
    cells = np.flatnonzero(values >= float(threshold))
    if cells.size == 0:
        return
    lon_min, lon_max, lat_min, lat_max = extent
    inside = (
        (grid.lon[cells] >= lon_min)
        & (grid.lon[cells] <= lon_max)
        & (grid.lat[cells] >= lat_min)
        & (grid.lat[cells] <= lat_max)
    )
    cells = cells[inside]
    if cells.size == 0:
        return
    if cells.size > int(max_count):
        cells = cells[np.argsort(values[cells])[-int(max_count):]]
    sizes = float(size_base) + float(size_scale) * np.clip(values[cells], 0.0, 1.0)
    ax.scatter(
        grid.lon[cells],
        grid.lat[cells],
        s=sizes,
        c=color,
        marker=marker,
        edgecolors=edgecolor,
        linewidths=0.35,
        alpha=0.92,
        zorder=zorder,
    )


def _render_coastal_morphology(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    plain_r = render.to_raster_continuous(
        grid, coastal.coastal_plain_rank, width=width, height=height)
    cliff_r = render.to_raster_continuous(
        grid, coastal.coastal_cliff_rank, width=width, height=height)
    shore_r = render.to_raster_continuous(
        grid, coastal.shoreface_rank, width=width, height=height)
    barrier_r = render.to_raster_continuous(
        grid, coastal.barrier_lagoon_rank, width=width, height=height)
    estuary_r = render.to_raster_continuous(
        grid, coastal.estuary_rank, width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 coastal morphology ranks (topology preserved)",
    )
    shore = np.ma.masked_where(shore_r <= 0.10, shore_r)
    ax.imshow(
        shore,
        cmap="PuBuGn",
        vmin=0.10,
        vmax=1.0,
        alpha=np.clip(shore_r * 0.38, 0.0, 0.42),
        extent=[-180, 180, -90, 90],
        zorder=3,
    )
    plain = np.ma.masked_where(plain_r <= 0.10, plain_r)
    ax.imshow(
        plain,
        cmap="Greens",
        vmin=0.10,
        vmax=1.0,
        alpha=np.clip(plain_r * 0.48, 0.0, 0.52),
        extent=[-180, 180, -90, 90],
        zorder=4,
    )
    estuary = np.ma.masked_where(estuary_r <= 0.10, estuary_r)
    ax.imshow(
        estuary,
        cmap="Blues",
        vmin=0.10,
        vmax=1.0,
        alpha=np.clip(estuary_r * 0.66, 0.0, 0.76),
        extent=[-180, 180, -90, 90],
        zorder=5,
    )
    barrier = np.ma.masked_where(barrier_r <= 0.10, barrier_r)
    ax.imshow(
        barrier,
        cmap="Wistia",
        vmin=0.10,
        vmax=1.0,
        alpha=np.clip(barrier_r * 0.66, 0.0, 0.78),
        extent=[-180, 180, -90, 90],
        zorder=6,
    )
    cliff = np.ma.masked_where(cliff_r <= 0.10, cliff_r)
    ax.imshow(
        cliff,
        cmap="Reds",
        vmin=0.10,
        vmax=1.0,
        alpha=np.clip(cliff_r * 0.62, 0.0, 0.74),
        extent=[-180, 180, -90, 90],
        zorder=7,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    plain_r = render.to_raster_continuous(
        grid, coastal.coastal_plain_rank, width=width, height=height)
    cliff_r = render.to_raster_continuous(
        grid, coastal.coastal_cliff_rank, width=width, height=height)
    shore_r = render.to_raster_continuous(
        grid, coastal.shoreface_rank, width=width, height=height)
    barrier_r = render.to_raster_continuous(
        grid, coastal.barrier_lagoon_rank, width=width, height=height)
    estuary_r = render.to_raster_continuous(
        grid, coastal.estuary_rank, width=width, height=height)
    centers = _coastal_zoom_centers(grid, coastal)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.55 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=38.0, lat_height=22.0)
        plain_crop, _ = _crop_lon_lat(
            plain_r, lon, lat, lon_width=38.0, lat_height=22.0)
        cliff_crop, _ = _crop_lon_lat(
            cliff_r, lon, lat, lon_width=38.0, lat_height=22.0)
        shore_crop, _ = _crop_lon_lat(
            shore_r, lon, lat, lon_width=38.0, lat_height=22.0)
        barrier_crop, _ = _crop_lon_lat(
            barrier_r, lon, lat, lon_width=38.0, lat_height=22.0)
        estuary_crop, _ = _crop_lon_lat(
            estuary_r, lon, lat, lon_width=38.0, lat_height=22.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        for crop, cmap, threshold, alpha_scale, vmax in (
            (shore_crop, "PuBuGn", 0.10, 0.52, 1.0),
            (plain_crop, "Greens", 0.10, 0.62, 1.0),
            (estuary_crop, "Blues", 0.10, 0.78, 1.0),
            (barrier_crop, "Wistia", 0.10, 0.78, 1.0),
            (cliff_crop, "Reds", 0.10, 0.76, 1.0),
        ):
            masked = np.ma.masked_where(crop <= threshold, crop)
            axes[row, 1].imshow(
                masked,
                cmap=cmap,
                vmin=threshold,
                vmax=vmax,
                alpha=np.clip(crop * alpha_scale, 0.0, min(0.92, alpha_scale + 0.12)),
                extent=extent,
            )
        axes[row, 1].set_title(f"{label}: plain/cliff/shore/barrier/estuary")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_process_linework(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    micro: _Microgeomorphology | None,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 coastal process linework",
    )
    _draw_coastal_process_linework(
        ax,
        grid,
        refined_rel,
        coastal,
        micro=micro,
        zoom=False,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_process_linework_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    micro: _Microgeomorphology | None,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    centers = _coastal_process_linework_zoom_centers(grid, coastal)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=32.0, lat_height=18.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: raster elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        _draw_coastal_process_linework(
            axes[row, 1],
            grid,
            refined_rel,
            coastal,
            micro=micro,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: distributary / estuary / barrier")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_process_microrelief_delta(
    grid: SphereGrid,
    coastal: _CoastalMorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        coastal.coastal_process_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 coastal process microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_process_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    micro: _Microgeomorphology | None,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        coastal.coastal_process_microrelief_delta,
        width=width,
        height=height,
    )
    centers = _coastal_process_linework_zoom_centers(grid, coastal)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        12.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=32.0, lat_height=18.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=32.0, lat_height=18.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        _draw_coastal_process_linework(
            axes[row, 1],
            grid,
            refined_rel,
            coastal,
            micro=micro,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: applied coastal process delta")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_depositional_microrelief_delta(
    grid: SphereGrid,
    coastal: _CoastalMorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    delta_r = render.to_raster_continuous(
        grid,
        coastal.coastal_depositional_microrelief_delta,
        width=width,
        height=height,
    )
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        10.0,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        delta_r,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title("72000 coastal depositional microrelief delta (m)")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coastal_depositional_microrelief_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid,
        coastal.coastal_depositional_microrelief_delta,
        width=width,
        height=height,
    )
    plain_r = render.to_raster_continuous(
        grid, coastal.coastal_depositional_plain_rank, width=width, height=height)
    strand_r = render.to_raster_continuous(
        grid, coastal.strandplain_rank, width=width, height=height)
    tidal_r = render.to_raster_continuous(
        grid, coastal.tidal_flat_rank, width=width, height=height)
    centers = _coastal_depositional_zoom_centers(grid, coastal)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    abs_delta = np.abs(delta_r[np.isfinite(delta_r)])
    abs_nonzero = abs_delta[abs_delta > 1.0e-6]
    vmax = max(
        float(np.nanpercentile(abs_nonzero, 95.0)) if abs_nonzero.size else 0.0,
        10.0,
    )
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=34.0, lat_height=20.0)
        delta_crop, _ = _crop_lon_lat(
            delta_r, lon, lat, lon_width=34.0, lat_height=20.0)
        plain_crop, _ = _crop_lon_lat(
            plain_r, lon, lat, lon_width=34.0, lat_height=20.0)
        strand_crop, _ = _crop_lon_lat(
            strand_r, lon, lat, lon_width=34.0, lat_height=20.0)
        tidal_crop, _ = _crop_lon_lat(
            tidal_r, lon, lat, lon_width=34.0, lat_height=20.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: refined elevation")
        im1 = axes[row, 1].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
        )
        for crop, cmap, threshold, alpha_scale in (
            (plain_crop, "Greens", 0.10, 0.36),
            (tidal_crop, "PuBu", 0.10, 0.48),
            (strand_crop, "YlOrBr", 0.10, 0.56),
        ):
            masked = np.ma.masked_where(crop <= threshold, crop)
            axes[row, 1].imshow(
                masked,
                cmap=cmap,
                vmin=threshold,
                vmax=1.0,
                alpha=np.clip(crop * alpha_scale, 0.0, min(0.90, alpha_scale + 0.12)),
                extent=extent,
            )
        axes[row, 1].set_title(f"{label}: depositional delta + ranks")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _draw_coastal_process_linework(
    ax,
    grid: SphereGrid,
    refined_rel: np.ndarray,
    coastal: _CoastalMorphology,
    *,
    micro: _Microgeomorphology | None = None,
    extent: list[float] | None = None,
    zoom: bool,
) -> None:
    land = np.asarray(refined_rel, dtype=np.float64) >= 0.0
    if micro is not None:
        _draw_coastal_process_mouth_paths(
            ax,
            grid,
            land,
            coastal,
            micro,
            extent=extent,
            zoom=zoom,
        )
    barrier_cells = _coastal_linework_cells(
        grid,
        coastal.barrier_spit_rank,
        min_score=0.34,
        max_count=92 if zoom else 46,
        extent=extent,
        spacing_passes=1 if zoom else 2,
    )
    estuary_cells = _coastal_linework_cells(
        grid,
        coastal.estuary_funnel_rank,
        min_score=0.46,
        max_count=64 if zoom else 34,
        extent=extent,
        spacing_passes=2,
    )
    delta_cells = _coastal_linework_cells(
        grid,
        coastal.delta_distributary_rank,
        min_score=0.52,
        max_count=72 if zoom else 38,
        extent=extent,
        spacing_passes=2,
    )
    mouth_bar_cells = _coastal_linework_cells(
        grid,
        coastal.delta_mouth_bar_rank,
        min_score=0.46,
        max_count=42 if zoom else 22,
        extent=extent,
        spacing_passes=2,
    )
    tidal_channel_cells = _coastal_linework_cells(
        grid,
        coastal.estuary_tidal_channel_rank,
        min_score=0.48,
        max_count=34 if zoom else 18,
        extent=extent,
        spacing_passes=2,
    )
    for cell in barrier_cells:
        _draw_barrier_spit_symbol(
            ax,
            grid,
            land,
            coastal,
            int(cell),
            zoom=zoom,
        )
    for cell in estuary_cells:
        _draw_estuary_funnel_symbol(
            ax,
            grid,
            land,
            coastal,
            int(cell),
            zoom=zoom,
        )
    for cell in tidal_channel_cells:
        _draw_estuary_tidal_channel_symbol(
            ax,
            grid,
            land,
            coastal,
            int(cell),
            zoom=zoom,
        )
    for cell in delta_cells:
        _draw_delta_distributary_symbol(
            ax,
            grid,
            land,
            coastal,
            int(cell),
            zoom=zoom,
        )
    for cell in mouth_bar_cells:
        _draw_delta_mouth_bar_symbol(
            ax,
            grid,
            land,
            coastal,
            int(cell),
            zoom=zoom,
        )


def _draw_coastal_process_mouth_paths(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    micro: _Microgeomorphology,
    *,
    extent: list[float] | None,
    zoom: bool,
) -> None:
    delta_rank = np.asarray(coastal.delta_distributary_rank, dtype=np.float64)
    estuary_rank = np.asarray(coastal.estuary_funnel_rank, dtype=np.float64)
    barrier_rank = np.asarray(coastal.barrier_spit_rank, dtype=np.float64)
    ocean_process = (delta_rank > 0.14) | (estuary_rank > 0.12)
    if not np.any(ocean_process):
        return
    river_rank = np.maximum.reduce(
        [
            np.asarray(micro.river_path_rank, dtype=np.float64),
            0.72 * np.asarray(micro.basin_trunk_rank, dtype=np.float64),
            0.42 * np.asarray(micro.river_rank, dtype=np.float64),
        ]
    )
    receiver = np.asarray(micro.river_receiver, dtype=np.int64)
    mouth_land = _land_neighbors_of(grid, ocean_process, land)
    mouth_land &= river_rank >= (0.14 if zoom else 0.20)
    mouths = np.flatnonzero(mouth_land)
    if mouths.size == 0:
        return
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        inside = (
            (grid.lon[mouths] >= lon_min)
            & (grid.lon[mouths] <= lon_max)
            & (grid.lat[mouths] >= lat_min)
            & (grid.lat[mouths] <= lat_max)
        )
        mouths = mouths[inside]
    if mouths.size == 0:
        return

    mouth_delta = np.zeros(grid.n, dtype=np.float64)
    mouth_estuary = np.zeros(grid.n, dtype=np.float64)
    mouth_barrier = np.zeros(grid.n, dtype=np.float64)
    for mouth in mouths:
        nbs = grid.neighbors[int(mouth)]
        ocean_nbs = nbs[~land[nbs]]
        if ocean_nbs.size == 0:
            continue
        mouth_delta[int(mouth)] = float(np.nanmax(delta_rank[ocean_nbs]))
        mouth_estuary[int(mouth)] = float(np.nanmax(estuary_rank[ocean_nbs]))
        mouth_barrier[int(mouth)] = float(np.nanmax(barrier_rank[ocean_nbs]))
    accumulation = np.asarray(micro.flow_accumulation, dtype=np.float64)
    acc = accumulation[mouths]
    acc_scale = max(float(np.nanpercentile(acc, 95.0)) if acc.size else 0.0, 1.0e-9)
    mouth_score = (
        0.50 * mouth_delta[mouths]
        + 0.28 * mouth_estuary[mouths]
        + 0.16 * river_rank[mouths]
        + 0.06 * np.clip(acc / acc_scale, 0.0, 1.0)
    )
    ordered = mouths[np.argsort(mouth_score)[::-1]]
    max_mouths = 5 if zoom else 30
    blocked = np.zeros(grid.n, dtype=bool)
    selected: list[int] = []
    for mouth in ordered:
        mouth = int(mouth)
        if blocked[mouth]:
            continue
        river_mouth_support = max(mouth_delta[mouth], mouth_estuary[mouth])
        if river_mouth_support <= (0.22 if zoom else 0.18):
            continue
        if mouth_barrier[mouth] > 0.16 and mouth_barrier[mouth] > 0.78 * river_mouth_support:
            continue
        selected.append(mouth)
        if len(selected) >= max_mouths:
            break
        _block_local_neighborhood(grid, blocked, mouth, passes=5 if zoom else 7)

    for mouth in selected:
        path = _coastal_mouth_upstream_path(
            grid,
            receiver,
            land,
            river_rank,
            accumulation,
            mouth,
            max_steps=6 if zoom else 5,
        )
        if len(path) >= 2:
            path = list(reversed(path))
            lons = grid.lon[np.asarray(path, dtype=np.int64)].astype(np.float64)
            lats = grid.lat[np.asarray(path, dtype=np.int64)].astype(np.float64)
            if not np.any(np.abs(np.diff(lons)) > 180.0):
                x, y = _smooth_river_path_xy(
                    int(path[0]),
                    lons,
                    lats,
                    float(np.nanmax(micro.meander_belt_rank[path])),
                    zoom=zoom,
                )
                strength = float(np.clip(np.nanmax(river_rank[path]), 0.0, 1.0))
                ax.plot(
                    x,
                    y,
                    color="#0568b0",
                    linewidth=(0.38 if zoom else 0.28) + (0.78 if zoom else 0.52) * strength,
                    alpha=0.34 + 0.26 * strength,
                    solid_capstyle="round",
                    solid_joinstyle="round",
                    zorder=8,
                )
        if mouth_delta[mouth] >= max(0.16, 0.82 * mouth_estuary[mouth]):
            _draw_delta_mouth_fan_symbol(
                ax,
                grid,
                land,
                coastal,
                mouth,
                zoom=zoom,
            )


def _coastal_mouth_upstream_path(
    grid: SphereGrid,
    receiver: np.ndarray,
    land: np.ndarray,
    river_rank: np.ndarray,
    accumulation: np.ndarray,
    mouth: int,
    *,
    max_steps: int,
) -> list[int]:
    path = [int(mouth)]
    current = int(mouth)
    seen = {current}
    for _ in range(max(0, int(max_steps))):
        nbs = grid.neighbors[current]
        candidates = nbs[
            land[nbs]
            & (receiver[nbs] == current)
            & (river_rank[nbs] >= 0.10)
        ]
        candidates = np.asarray([int(c) for c in candidates if int(c) not in seen], dtype=np.int64)
        if candidates.size == 0:
            break
        score = (
            0.68 * np.asarray(river_rank[candidates], dtype=np.float64)
            + 0.32
            * np.clip(
                np.asarray(accumulation[candidates], dtype=np.float64)
                / max(float(np.nanmax(accumulation[candidates])), 1.0e-9),
                0.0,
                1.0,
            )
        )
        nxt = int(candidates[int(np.nanargmax(score))])
        path.append(nxt)
        seen.add(nxt)
        current = nxt
    return path


def _draw_delta_mouth_fan_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    mouth: int,
    *,
    zoom: bool,
) -> None:
    nbs = grid.neighbors[int(mouth)]
    ocean_nbs = nbs[~land[nbs]]
    if ocean_nbs.size == 0:
        return
    delta_rank = np.asarray(coastal.delta_distributary_rank, dtype=np.float64)
    delta_cell = int(ocean_nbs[int(np.nanargmax(delta_rank[ocean_nbs]))])
    rank = float(np.clip(delta_rank[delta_cell], 0.0, 1.0))
    if rank <= 0.0:
        return
    lon = float(grid.lon[delta_cell])
    lat = float(grid.lat[delta_cell])
    angle = _coastal_normal_angle(grid, land, delta_cell)
    length = (0.98 if zoom else 1.24) * (0.72 + 0.54 * rank)
    for branch in (-24.0, 0.0, 24.0):
        x, y = _coastal_segment_xy(
            lon,
            lat,
            angle + branch,
            length * (0.76 if branch else 1.05),
            bend=0.10 * np.sign(branch),
        )
        ax.plot(
            x,
            y,
            color="#17c2ff" if branch else "#049ee4",
            linewidth=((0.72 if zoom else 0.48) + 1.25 * rank) * (0.78 if branch else 1.0),
            alpha=0.50 + 0.38 * rank,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=11,
        )


def _coastal_linework_cells(
    grid: SphereGrid,
    rank: np.ndarray,
    *,
    min_score: float,
    max_count: int,
    extent: list[float] | None,
    spacing_passes: int,
) -> np.ndarray:
    values = np.asarray(rank, dtype=np.float64)
    cells = np.flatnonzero(values >= float(min_score))
    if cells.size == 0:
        return cells
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        inside = (
            (grid.lon[cells] >= lon_min)
            & (grid.lon[cells] <= lon_max)
            & (grid.lat[cells] >= lat_min)
            & (grid.lat[cells] <= lat_max)
        )
        cells = cells[inside]
    if cells.size == 0:
        return cells
    ordered = cells[np.argsort(values[cells])[::-1]]
    selected: list[int] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for cell in ordered:
        cell = int(cell)
        if blocked[cell]:
            continue
        selected.append(cell)
        if len(selected) >= int(max_count):
            break
        if int(spacing_passes) > 0:
            seed = np.zeros(grid.n, dtype=bool)
            seed[cell] = True
            blocked |= _dilate_mask(grid, seed, passes=int(spacing_passes))
    return np.asarray(selected, dtype=np.int64)


def _draw_delta_distributary_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(coastal.delta_distributary_rank[cell], 0.0, 1.0))
    angle = _coastal_normal_angle(grid, land, cell)
    length = (0.72 if zoom else 0.95) * (0.72 + 0.58 * rank)
    linewidth = (0.72 if zoom else 0.52) + 1.05 * rank
    branches = (-18.0, 0.0, 18.0) if rank >= 0.48 else (-14.0, 14.0)
    for branch in branches:
        x, y = _coastal_segment_xy(
            lon,
            lat,
            angle + branch,
            length * (0.78 if branch else 1.0),
            bend=0.08 * np.sign(branch),
        )
        ax.plot(
            x,
            y,
            color="#37b7ff",
            linewidth=linewidth * (0.82 if branch else 1.0),
            alpha=0.48 + 0.42 * rank,
            solid_capstyle="round",
            zorder=10,
        )


def _draw_delta_mouth_bar_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(coastal.delta_mouth_bar_rank[cell], 0.0, 1.0))
    if rank <= 0.0:
        return
    normal = _coastal_normal_angle(grid, land, cell)
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    width = (0.46 if zoom else 0.62) * lon_scale * (0.72 + 0.42 * rank)
    height = (0.12 if zoom else 0.16) * (0.70 + 0.34 * rank)
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        width,
        height,
        angle=normal + 90.0,
        facecolor="#ffd17a",
        edgecolor="#d8861d",
        linewidth=0.60 if zoom else 0.42,
        alpha=0.40 + 0.34 * rank,
        zorder=12,
    ))


def _draw_estuary_tidal_channel_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(coastal.estuary_tidal_channel_rank[cell], 0.0, 1.0))
    if rank <= 0.0:
        return
    angle = _coastal_normal_angle(grid, land, cell)
    length = (0.58 if zoom else 0.78) * (0.70 + 0.48 * rank)
    x, y = _coastal_segment_xy(
        lon,
        lat,
        angle,
        length,
        bend=0.05,
    )
    ax.plot(
        x,
        y,
        color="#073b73",
        linewidth=(0.48 if zoom else 0.34) + 0.86 * rank,
        alpha=0.44 + 0.36 * rank,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=12,
    )


def _draw_estuary_funnel_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(coastal.estuary_funnel_rank[cell], 0.0, 1.0))
    angle = _coastal_normal_angle(grid, land, cell)
    length = (0.82 if zoom else 1.05) * (0.68 + 0.48 * rank)
    spread = 20.0 + 18.0 * rank
    for offset in (-spread, spread):
        x, y = _coastal_segment_xy(
            lon,
            lat,
            angle + offset,
            length,
            bend=0.05 * np.sign(offset),
        )
        ax.plot(
            x,
            y,
            color="#1d6ea8",
            linewidth=0.55 + 1.25 * rank,
            alpha=0.42 + 0.36 * rank,
            solid_capstyle="round",
            zorder=9,
        )
    x, y = _coastal_segment_xy(lon, lat, angle, length * 0.92, bend=0.0)
    ax.plot(
        x,
        y,
        color="#6fd8ff",
        linewidth=0.52 + 0.95 * rank,
        alpha=0.42 + 0.34 * rank,
        solid_capstyle="round",
        zorder=10,
    )


def _draw_barrier_spit_symbol(
    ax,
    grid: SphereGrid,
    land: np.ndarray,
    coastal: _CoastalMorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(coastal.barrier_spit_rank[cell], 0.0, 1.0))
    normal = _coastal_normal_angle(grid, land, cell)
    angle = normal + 90.0
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    width = (0.62 if zoom else 0.82) * lon_scale * (0.86 + 0.52 * rank)
    height = (0.13 if zoom else 0.18) * (0.76 + 0.42 * rank)
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        width,
        height,
        angle=angle,
        facecolor="#f3df88",
        edgecolor="#8d7939",
        linewidth=0.72 if zoom else 0.48,
        alpha=0.62 + 0.26 * rank,
        zorder=8,
    ))


def _coastal_normal_angle(
    grid: SphereGrid,
    land: np.ndarray,
    cell: int,
) -> float:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    nbs = grid.neighbors[int(cell)]
    land_nbs = nbs[np.asarray(land, dtype=bool)[nbs]]
    if land_nbs.size:
        lon_values = np.asarray(grid.lon[land_nbs], dtype=np.float64)
        lon_values = lon + ((lon_values - lon + 180.0) % 360.0 - 180.0)
        dx = lon - float(np.nanmean(lon_values))
        dy = lat - float(np.nanmean(grid.lat[land_nbs]))
        if abs(dx) + abs(dy) > 1.0e-9:
            return float(np.degrees(np.arctan2(dy, dx)))
    return _microshape_angle(cell)


def _coastal_segment_xy(
    lon: float,
    lat: float,
    angle_deg: float,
    length: float,
    *,
    bend: float,
) -> tuple[np.ndarray, np.ndarray]:
    theta = np.radians(angle_deg)
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    dx = float(length) * lon_scale * np.cos(theta)
    dy = float(length) * np.sin(theta)
    cross_x = -np.sin(theta) * float(length) * lon_scale * float(bend)
    cross_y = np.cos(theta) * float(length) * float(bend)
    t = np.linspace(0.0, 1.0, 9)
    curve = np.sin(np.pi * t)
    x = lon + dx * t + cross_x * curve
    y = lat + dy * t + cross_y * curve
    return x, y


def _render_hydrology(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    lake_r = render.to_raster(
        grid, micro.lake_mask.astype(np.float64), width=width, height=height)
    lake_basin_r = render.to_raster_continuous(
        grid, micro.lake_basin_rank, width=width, height=height)
    delta_fan_r = render.to_raster_continuous(
        grid, micro.delta_fan_rank, width=width, height=height)
    delta_plain_r = render.to_raster(
        grid, micro.delta_plain_mask.astype(np.float64), width=width, height=height)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(ax, elev_r, title="72000 selected-snapshot hydrology")
    lake = np.ma.masked_where(lake_r <= 0.5, lake_r)
    lake_basin = np.ma.masked_where(lake_basin_r <= 0.08, lake_basin_r)
    ax.imshow(
        lake_basin,
        cmap="Blues",
        vmin=0.08,
        vmax=1.0,
        alpha=np.clip(lake_basin_r * 0.40, 0.0, 0.45),
        extent=[-180, 180, -90, 90],
    )
    ax.imshow(
        lake,
        cmap="winter",
        vmin=0.0,
        vmax=1.0,
        alpha=0.85,
        extent=[-180, 180, -90, 90],
    )
    delta = np.ma.masked_where(np.maximum(delta_fan_r, delta_plain_r) <= 0.05,
                               np.maximum(delta_fan_r, delta_plain_r))
    ax.imshow(
        delta,
        cmap="autumn",
        vmin=0.0,
        vmax=1.0,
        alpha=np.clip(np.maximum(delta_fan_r, delta_plain_r) * 0.68, 0.0, 0.72),
        extent=[-180, 180, -90, 90],
    )
    _draw_river_segments(
        ax,
        grid,
        micro,
        min_rank=0.94,
        rank_values=micro.river_rank,
        color="#2a82bd",
        exclude_mask=micro.river_path_rank > 0.0,
        linewidth_base=0.18,
        linewidth_scale=0.60,
        alpha_base=0.16,
        alpha_scale=0.30,
        max_segments=140,
        spacing_passes=2,
    )
    _draw_river_paths(
        ax,
        grid,
        micro,
        min_rank=0.38,
        rank_values=micro.river_path_rank,
        meander_values=micro.meander_belt_rank,
        color="#053f8f",
        exclude_mask=micro.basin_trunk_rank > 0.0,
        linewidth_base=0.32,
        linewidth_scale=1.55,
        alpha_base=0.42,
        alpha_scale=0.52,
        min_path_cells=6,
        max_paths=48,
    )
    _draw_river_paths(
        ax,
        grid,
        micro,
        min_rank=0.18,
        rank_values=micro.basin_trunk_rank,
        meander_values=micro.meander_belt_rank,
        color="#043985",
        linewidth_base=0.42,
        linewidth_scale=1.62,
        alpha_base=0.52,
        alpha_scale=0.42,
        min_path_cells=4,
        max_paths=64,
    )
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_hydrology_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    lake_r = render.to_raster(
        grid, micro.lake_mask.astype(np.float64), width=width, height=height)
    lake_basin_r = render.to_raster_continuous(
        grid, micro.lake_basin_rank, width=width, height=height)
    delta_fan_r = render.to_raster_continuous(
        grid, micro.delta_fan_rank, width=width, height=height)
    delta_plain_r = render.to_raster(
        grid, micro.delta_plain_mask.astype(np.float64), width=width, height=height)
    hydro_r = np.maximum.reduce([
        lake_basin_r * 0.82,
        lake_r * 1.15,
        delta_fan_r * 1.05,
        delta_plain_r * 0.92,
    ])
    centers = _hydrology_zoom_centers(grid, refined_rel, micro)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.9 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(elev_r, lon, lat, lon_width=54.0, lat_height=28.0)
        hydro_crop, _ = _crop_lon_lat(hydro_r, lon, lat, lon_width=54.0, lat_height=28.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        hydro = np.ma.masked_where(hydro_crop <= 0.05, hydro_crop)
        axes[row, 1].imshow(
            hydro,
            cmap="Blues",
            vmin=0.05,
            vmax=1.15,
            alpha=np.clip(hydro_crop * 0.72, 0.0, 0.88),
            extent=extent,
        )
        _draw_lake_shoreline_linework(
            axes[row, 1],
            grid,
            micro,
            extent=extent,
            zoom=True,
        )
        segment_limit = 14 if label == "lake basin" else 24
        path_limit = 8 if label == "lake basin" else 18
        _draw_river_segments(
            axes[row, 1],
            grid,
            micro,
            min_rank=0.86,
            rank_values=micro.river_rank,
            color="#2a82bd",
            exclude_mask=micro.river_path_rank > 0.0,
            linewidth_base=0.35,
            linewidth_scale=1.25,
            alpha_base=0.24,
            alpha_scale=0.38,
            extent=extent,
            zoom=True,
            max_segments=segment_limit,
            spacing_passes=1,
        )
        _draw_river_paths(
            axes[row, 1],
            grid,
            micro,
            min_rank=0.34,
            rank_values=micro.river_path_rank,
            meander_values=micro.meander_belt_rank,
            color="#053f8f",
            exclude_mask=micro.basin_trunk_rank > 0.0,
            linewidth_base=0.75,
            linewidth_scale=2.55,
            alpha_base=0.52,
            alpha_scale=0.44,
            extent=extent,
            zoom=True,
            min_path_cells=4,
            max_paths=path_limit,
        )
        _draw_river_paths(
            axes[row, 1],
            grid,
            micro,
            min_rank=0.16,
            rank_values=micro.basin_trunk_rank,
            meander_values=micro.meander_belt_rank,
            color="#043985",
            linewidth_base=0.85,
            linewidth_scale=2.65,
            alpha_base=0.58,
            alpha_scale=0.36,
            extent=extent,
            zoom=True,
            min_path_cells=3,
            max_paths=path_limit,
        )
        axes[row, 1].set_title(f"{label}: rivers/lakes/deltas")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.75)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_fluvial_lacustrine_microshapes(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elev_r,
        title="72000 fluvial / lacustrine microshapes",
    )
    _draw_fluvial_lacustrine_microshapes(ax, grid, micro, zoom=False)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_fluvial_lacustrine_microshapes_zoom_sheet(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    elev_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    centers = _fluvial_lacustrine_zoom_centers(grid, micro)
    fig, axes = plt.subplots(
        len(centers),
        2,
        figsize=(10.5, 3.6 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    for row, (label, lon, lat) in enumerate(centers):
        elev_crop, extent = _crop_lon_lat(
            elev_r, lon, lat, lon_width=34.0, lat_height=20.0)
        im0 = axes[row, 0].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: raster elevation")
        axes[row, 1].imshow(
            elev_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        _draw_fluvial_lacustrine_microshapes(
            axes[row, 1],
            grid,
            micro,
            extent=extent,
            zoom=True,
        )
        axes[row, 1].set_title(f"{label}: scrolls / swales / shorelines")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.72)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _draw_fluvial_lacustrine_microshapes(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    *,
    extent: list[float] | None = None,
    zoom: bool,
) -> None:
    swale_cells = _hydrology_microshape_cells(
        grid,
        micro.floodplain_swale_rank,
        min_score=0.20,
        max_count=180 if zoom else 72,
        extent=extent,
        spacing_passes=1 if zoom else 2,
    )
    scroll_cells = _hydrology_microshape_cells(
        grid,
        micro.meander_scroll_rank,
        min_score=0.09,
        max_count=160 if zoom else 64,
        extent=extent,
        spacing_passes=2,
    )
    lake_cells = _hydrology_microshape_cells(
        grid,
        micro.lake_shoreline_rank,
        min_score=0.20,
        max_count=56 if zoom else 28,
        extent=extent,
        spacing_passes=2,
    )
    _draw_lake_shoreline_linework(ax, grid, micro, extent=extent, zoom=zoom)
    for cell in swale_cells:
        _draw_floodplain_swale_symbol(ax, grid, micro, int(cell), zoom=zoom)
    for cell in scroll_cells:
        _draw_meander_scroll_symbol(ax, grid, micro, int(cell), zoom=zoom)
    for cell in lake_cells:
        _draw_lake_shoreline_symbol(ax, grid, micro, int(cell), zoom=zoom)
    _draw_river_paths(
        ax,
        grid,
        micro,
        min_rank=0.30 if zoom else 0.34,
        rank_values=micro.basin_trunk_rank,
        meander_values=micro.meander_belt_rank,
        color="#064b9b",
        linewidth_base=0.36 if zoom else 0.30,
        linewidth_scale=1.22 if zoom else 0.92,
        alpha_base=0.44,
        alpha_scale=0.36,
        extent=extent,
        zoom=zoom,
        min_path_cells=4 if zoom else 3,
        max_paths=18 if zoom else 72,
    )


def _hydrology_microshape_cells(
    grid: SphereGrid,
    rank: np.ndarray,
    *,
    min_score: float,
    max_count: int,
    extent: list[float] | None,
    spacing_passes: int,
) -> np.ndarray:
    values = np.asarray(rank, dtype=np.float64)
    cells = np.flatnonzero(values >= float(min_score))
    if cells.size == 0:
        return cells
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        inside = (
            (grid.lon[cells] >= lon_min)
            & (grid.lon[cells] <= lon_max)
            & (grid.lat[cells] >= lat_min)
            & (grid.lat[cells] <= lat_max)
        )
        cells = cells[inside]
    if cells.size == 0:
        return cells
    ordered = cells[np.argsort(values[cells])[::-1]]
    selected: list[int] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for cell in ordered:
        cell = int(cell)
        if blocked[cell]:
            continue
        selected.append(cell)
        if len(selected) >= int(max_count):
            break
        if int(spacing_passes) > 0:
            seed = np.zeros(grid.n, dtype=bool)
            seed[cell] = True
            blocked |= _dilate_mask(grid, seed, passes=int(spacing_passes))
    return np.asarray(selected, dtype=np.int64)


def _draw_meander_scroll_symbol(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(micro.meander_scroll_rank[cell], 0.0, 1.0))
    angle = _river_local_angle(grid, micro, cell)
    length = (0.80 if zoom else 1.02) * (0.60 + 0.48 * rank)
    offset = (0.18 if zoom else 0.24) * (0.45 + 0.40 * rank)
    for side in (-1.0, 1.0):
        x, y = _hydrology_segment_xy(
            lon,
            lat,
            angle + 7.0 * side,
            length,
            lateral=offset * side,
            wiggle=0.18 * side,
        )
        ax.plot(
            x,
            y,
            color="#65a85b",
            linewidth=0.48 + 1.10 * rank,
            alpha=0.34 + 0.42 * rank,
            solid_capstyle="round",
            zorder=7,
        )


def _draw_floodplain_swale_symbol(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(micro.floodplain_swale_rank[cell], 0.0, 1.0))
    angle = _river_local_angle(grid, micro, cell) + 90.0
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    width = (0.42 if zoom else 0.58) * lon_scale * (0.74 + 0.46 * rank)
    height = (0.10 if zoom else 0.14) * (0.70 + 0.40 * rank)
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        width,
        height,
        angle=angle,
        facecolor="#87c978",
        edgecolor="#4d8b50",
        linewidth=0.42 if zoom else 0.28,
        alpha=0.24 + 0.32 * rank,
        zorder=6,
    ))


def _draw_lake_shoreline_symbol(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    cell: int,
    *,
    zoom: bool,
) -> None:
    lon = float(grid.lon[cell])
    lat = float(grid.lat[cell])
    rank = float(np.clip(micro.lake_shoreline_rank[cell], 0.0, 1.0))
    angle = _microshape_angle(cell)
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    width = (0.54 if zoom else 0.72) * lon_scale * (0.80 + 0.48 * rank)
    height = (0.28 if zoom else 0.36) * (0.70 + 0.38 * rank)
    ax.add_patch(matplotlib.patches.Ellipse(
        (lon, lat),
        width,
        height,
        angle=angle,
        facecolor=(0.25, 0.78, 0.92, 0.20 + 0.30 * rank),
        edgecolor="#096a94",
        linewidth=0.72 if zoom else 0.46,
        alpha=0.66 + 0.22 * rank,
        zorder=8,
    ))


def _draw_lake_shoreline_linework(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    *,
    extent: list[float] | None,
    zoom: bool,
) -> None:
    rank = np.asarray(micro.lake_shoreline_rank, dtype=np.float64)
    active = rank >= (0.14 if zoom else 0.18)
    lake_adjacent = np.zeros(grid.n, dtype=bool)
    for lake_cell in np.flatnonzero(np.asarray(micro.lake_mask, dtype=bool)):
        lake_adjacent[grid.neighbors[int(lake_cell)]] = True
    if np.any(active & lake_adjacent):
        active &= lake_adjacent
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        active &= (
            (grid.lon >= lon_min)
            & (grid.lon <= lon_max)
            & (grid.lat >= lat_min)
            & (grid.lat <= lat_max)
        )
    if not np.any(active):
        return
    seen = np.zeros(grid.n, dtype=bool)
    components: list[tuple[float, np.ndarray]] = []
    for start in np.flatnonzero(active):
        start = int(start)
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        cells: list[int] = []
        while stack:
            cell = stack.pop()
            cells.append(cell)
            for nb in grid.neighbors[cell]:
                nb = int(nb)
                if active[nb] and not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        if len(cells) < 2:
            continue
        comp = np.asarray(cells, dtype=np.int64)
        score = float(np.nanmax(rank[comp]) * np.sqrt(float(comp.size)))
        components.append((score, comp))
    if not components:
        return
    components.sort(key=lambda item: item[0], reverse=True)
    max_components = 10 if zoom else 28
    max_edges_per_component = 90 if zoom else 34
    for _, comp in components[:max_components]:
        comp_mask = np.zeros(grid.n, dtype=bool)
        comp_mask[comp] = True
        edges: list[tuple[float, int, int]] = []
        for cell in comp:
            cell = int(cell)
            lon0 = float(grid.lon[cell])
            for nb in grid.neighbors[cell]:
                nb = int(nb)
                if nb <= cell or not comp_mask[nb]:
                    continue
                lon1 = float(grid.lon[nb])
                if abs(lon0 - lon1) > 180.0:
                    continue
                edge_score = 0.5 * (float(rank[cell]) + float(rank[nb]))
                edges.append((edge_score, cell, nb))
        if not edges:
            continue
        if len(edges) > max_edges_per_component:
            edges.sort(key=lambda item: item[0], reverse=True)
            edges = edges[:max_edges_per_component]
        for edge_score, cell, nb in edges:
            strength = float(np.clip(edge_score, 0.0, 1.0))
            ax.plot(
                [float(grid.lon[cell]), float(grid.lon[nb])],
                [float(grid.lat[cell]), float(grid.lat[nb])],
                color="#064f75",
                linewidth=(0.56 if zoom else 0.36) + (1.04 if zoom else 0.58) * strength,
                alpha=0.34 + 0.38 * strength,
                solid_capstyle="round",
                zorder=7,
            )


def _river_local_angle(
    grid: SphereGrid,
    micro: _Microgeomorphology,
    cell: int,
) -> float:
    rec = int(micro.river_receiver[int(cell)])
    lon = float(grid.lon[int(cell)])
    lat = float(grid.lat[int(cell)])
    if rec >= 0 and rec != int(cell):
        lon1 = float(grid.lon[rec])
        if abs(lon1 - lon) <= 180.0:
            dx = lon1 - lon
            dy = float(grid.lat[rec]) - lat
            if abs(dx) + abs(dy) > 1.0e-9:
                return float(np.degrees(np.arctan2(dy, dx)))
    nbs = grid.neighbors[int(cell)]
    river_nbs = nbs[np.asarray(micro.basin_trunk_rank[nbs] > 0.0)]
    if river_nbs.size:
        lon_values = np.asarray(grid.lon[river_nbs], dtype=np.float64)
        lon_values = lon + ((lon_values - lon + 180.0) % 360.0 - 180.0)
        dx = float(np.nanmean(lon_values)) - lon
        dy = float(np.nanmean(grid.lat[river_nbs])) - lat
        if abs(dx) + abs(dy) > 1.0e-9:
            return float(np.degrees(np.arctan2(dy, dx)))
    return _microshape_angle(cell)


def _hydrology_segment_xy(
    lon: float,
    lat: float,
    angle_deg: float,
    length: float,
    *,
    lateral: float,
    wiggle: float,
) -> tuple[np.ndarray, np.ndarray]:
    theta = np.radians(angle_deg)
    lon_scale = 1.0 / max(float(np.cos(np.radians(lat))), 0.45)
    dx = float(length) * lon_scale * np.cos(theta)
    dy = float(length) * np.sin(theta)
    normal_x = -np.sin(theta) * lon_scale
    normal_y = np.cos(theta)
    t = np.linspace(-0.5, 0.5, 13)
    curve = np.sin(2.0 * np.pi * (t + 0.5))
    x = lon + dx * t + normal_x * float(lateral) + normal_x * float(wiggle) * curve
    y = lat + dy * t + normal_y * float(lateral) + normal_y * float(wiggle) * curve
    return x, y


def _draw_river_segments(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    *,
    min_rank: float,
    rank_values: np.ndarray | None = None,
    meander_values: np.ndarray | None = None,
    color: str = "#0757a6",
    exclude_mask: np.ndarray | None = None,
    linewidth_base: float = 0.25,
    linewidth_scale: float = 1.25,
    alpha_base: float = 0.30,
    alpha_scale: float = 0.55,
    extent: list[float] | None = None,
    zoom: bool = False,
    max_segments: int | None = None,
    spacing_passes: int = 0,
) -> None:
    if rank_values is None:
        rank_values = micro.river_rank
    rank = np.asarray(rank_values, dtype=np.float64)
    meander = None
    if meander_values is not None:
        meander = np.asarray(meander_values, dtype=np.float64)
    receiver = np.asarray(micro.river_receiver, dtype=np.int64)
    cells = np.flatnonzero((rank >= float(min_rank)) & (receiver >= 0))
    if exclude_mask is not None and cells.size:
        excluded = np.asarray(exclude_mask, dtype=bool)
        cells = cells[~excluded[cells]]
    if extent is not None:
        lon_min, lon_max, lat_min, lat_max = extent
        inside = (
            (grid.lon[cells] >= lon_min)
            & (grid.lon[cells] <= lon_max)
            & (grid.lat[cells] >= lat_min)
            & (grid.lat[cells] <= lat_max)
        )
        cells = cells[inside]
    cells = cells[np.argsort(rank[cells])[::-1]]
    if max_segments is not None:
        selected: list[int] = []
        blocked = np.zeros(rank.shape, dtype=bool)
        for cell in cells:
            cell = int(cell)
            if blocked[cell]:
                continue
            selected.append(cell)
            if len(selected) >= int(max_segments):
                break
            if spacing_passes > 0:
                _block_local_neighborhood(
                    grid,
                    blocked,
                    cell,
                    passes=int(spacing_passes),
                )
        cells = np.asarray(selected, dtype=np.int64)
    for cell in cells:
        rec = int(receiver[int(cell)])
        if rec < 0 or rec == int(cell):
            continue
        lon0 = float(grid.lon[int(cell)])
        lon1 = float(grid.lon[rec])
        if abs(lon0 - lon1) > 180.0:
            continue
        lat0 = float(grid.lat[int(cell)])
        lat1 = float(grid.lat[rec])
        strength = float(np.clip(rank[int(cell)], 0.0, 1.0))
        linewidth = linewidth_base + linewidth_scale * strength
        alpha = alpha_base + alpha_scale * strength
        if zoom:
            linewidth *= 1.18
            alpha = min(0.98, alpha + 0.05)
        x, y = _river_segment_xy(
            int(cell),
            lon0,
            lat0,
            lon1,
            lat1,
            float(meander[int(cell)]) if meander is not None else 0.0,
            zoom=zoom,
        )
        ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            solid_capstyle="round",
            zorder=5,
        )


def _draw_river_paths(
    ax,
    grid: SphereGrid,
    micro: _Microgeomorphology,
    *,
    min_rank: float,
    rank_values: np.ndarray,
    meander_values: np.ndarray | None = None,
    color: str,
    linewidth_base: float,
    linewidth_scale: float,
    alpha_base: float,
    alpha_scale: float,
    extent: list[float] | None = None,
    zoom: bool = False,
    exclude_mask: np.ndarray | None = None,
    min_path_cells: int = 2,
    max_paths: int | None = None,
) -> None:
    rank = np.asarray(rank_values, dtype=np.float64)
    receiver = np.asarray(micro.river_receiver, dtype=np.int64)
    meander = (
        np.asarray(meander_values, dtype=np.float64)
        if meander_values is not None
        else np.zeros_like(rank)
    )
    active = rank >= float(min_rank)
    if exclude_mask is not None:
        active &= ~np.asarray(exclude_mask, dtype=bool)
    if not np.any(active):
        return
    active_cells = np.flatnonzero(active)
    has_upstream = np.zeros(rank.shape, dtype=bool)
    for cell in active_cells:
        rec = int(receiver[int(cell)])
        if rec >= 0 and rec != int(cell) and active[rec]:
            has_upstream[rec] = True
    sources = active_cells[~has_upstream[active_cells]]
    if sources.size == 0:
        sources = active_cells[np.argsort(rank[active_cells])[::-1][:64]]
    sources = sources[np.argsort(rank[sources])[::-1]]
    drawn = np.zeros(rank.shape, dtype=bool)
    path_limit = int(max_paths) if max_paths is not None else (160 if zoom else 96)
    for source in sources[:path_limit]:
        path: list[int] = []
        cell = int(source)
        seen: set[int] = set()
        for _ in range(512):
            if cell in seen or cell < 0 or not active[cell]:
                break
            path.append(cell)
            seen.add(cell)
            rec = int(receiver[cell])
            if rec < 0 or rec == cell:
                break
            cell = rec
        if len(path) < int(min_path_cells):
            continue
        if np.count_nonzero(drawn[np.asarray(path, dtype=np.int64)]) > max(1, len(path) // 3):
            continue
        lons = grid.lon[np.asarray(path, dtype=np.int64)].astype(np.float64)
        lats = grid.lat[np.asarray(path, dtype=np.int64)].astype(np.float64)
        if np.any(np.abs(np.diff(lons)) > 180.0):
            continue
        if extent is not None:
            lon_min, lon_max, lat_min, lat_max = extent
            if not (
                np.any((lons >= lon_min) & (lons <= lon_max))
                and np.any((lats >= lat_min) & (lats <= lat_max))
            ):
                continue
        strength = float(np.clip(np.nanmax(rank[path]), 0.0, 1.0))
        meander_strength = float(np.clip(np.nanmax(meander[path]), 0.0, 1.0))
        x, y = _smooth_river_path_xy(
            int(source),
            lons,
            lats,
            meander_strength,
            zoom=zoom,
        )
        linewidth = linewidth_base + linewidth_scale * strength
        alpha = alpha_base + alpha_scale * strength
        if zoom:
            linewidth *= 1.10
            alpha = min(0.98, alpha + 0.05)
        ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=6,
        )
        drawn[np.asarray(path, dtype=np.int64)] = True


def _block_local_neighborhood(
    grid: SphereGrid,
    blocked: np.ndarray,
    cell: int,
    *,
    passes: int,
) -> None:
    blocked[int(cell)] = True
    frontier = np.asarray([int(cell)], dtype=np.int64)
    for _ in range(max(0, int(passes))):
        if frontier.size == 0:
            break
        nbs = np.unique(np.concatenate([grid.neighbors[int(c)] for c in frontier]))
        nbs = nbs[~blocked[nbs]]
        if nbs.size == 0:
            break
        blocked[nbs] = True
        frontier = nbs.astype(np.int64, copy=False)


def _smooth_river_path_xy(
    source: int,
    lons: np.ndarray,
    lats: np.ndarray,
    meander_rank: float,
    *,
    zoom: bool,
) -> tuple[np.ndarray, np.ndarray]:
    pts = np.column_stack([lons, lats]).astype(np.float64)
    if pts.shape[0] <= 2:
        return _river_segment_xy(
            int(source),
            float(pts[0, 0]),
            float(pts[0, 1]),
            float(pts[-1, 0]),
            float(pts[-1, 1]),
            meander_rank,
            zoom=zoom,
        )
    smooth = pts.copy()
    iterations = 2 if zoom else 1
    for _ in range(iterations):
        if smooth.shape[0] <= 2:
            break
        left = smooth[:-1] * 0.75 + smooth[1:] * 0.25
        right = smooth[:-1] * 0.25 + smooth[1:] * 0.75
        new = [smooth[0]]
        for a, b in zip(left, right, strict=False):
            new.append(a)
            new.append(b)
        new.append(smooth[-1])
        smooth = np.asarray(new, dtype=np.float64)
    strength = float(np.sqrt(np.clip(meander_rank, 0.0, 1.0)))
    if smooth.shape[0] >= (8 if zoom else 6):
        # Long rendered paths should not read as rigid receiver-tree polylines
        # even where the terminal meander rank is low.
        strength = max(strength, 0.13 if zoom else 0.08)
    if strength > 0.05 and smooth.shape[0] >= 4:
        phase = ((int(source) * 1664525 + 1013904223) & 0xFFFF) / 65535.0
        diffs = np.gradient(smooth, axis=0)
        lengths = np.maximum(np.linalg.norm(diffs, axis=1), 1.0e-9)
        normals = np.column_stack([-diffs[:, 1] / lengths, diffs[:, 0] / lengths])
        t = np.linspace(0.0, 1.0, smooth.shape[0])
        envelope = np.sin(np.pi * t)
        wave = np.sin((2.0 + phase) * np.pi * t + 2.0 * np.pi * phase)
        mean_step = float(np.nanmedian(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        amp = min(0.80 if zoom else 0.46, max(0.05, 0.22 * mean_step)) * strength
        smooth = smooth + normals * (amp * envelope * wave)[:, None]
    return smooth[:, 0], smooth[:, 1]


def _river_segment_xy(
    cell: int,
    lon0: float,
    lat0: float,
    lon1: float,
    lat1: float,
    meander_rank: float,
    *,
    zoom: bool,
) -> tuple[np.ndarray, np.ndarray]:
    strength = float(np.sqrt(np.clip(meander_rank, 0.0, 1.0)))
    if strength <= 0.05:
        return np.asarray([lon0, lon1]), np.asarray([lat0, lat1])
    dx = lon1 - lon0
    dy = lat1 - lat0
    length = float(np.hypot(dx, dy))
    if length <= 1.0e-9:
        return np.asarray([lon0, lon1]), np.asarray([lat0, lat1])
    nx = -dy / length
    ny = dx / length
    # Deterministic but local phase; endpoints remain pinned to receiver cells.
    phase = ((int(cell) * 1103515245 + 12345) & 0xFFFF) / 65535.0
    side = -1.0 if phase < 0.5 else 1.0
    t = np.linspace(0.0, 1.0, 7 if zoom else 5)
    envelope = np.sin(np.pi * t)
    secondary = 0.42 * np.sin(2.0 * np.pi * t + 2.0 * np.pi * phase)
    amp = min(0.72 if zoom else 0.42, max(0.06, 0.30 * length)) * strength
    offset = side * amp * envelope * (1.0 + secondary)
    x = lon0 + dx * t + nx * offset
    y = lat0 + dy * t + ny * offset
    return x, y


def _hydrology_zoom_centers(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
) -> list[tuple[str, float, float]]:
    lat_weight = _zoom_temperate_latitude_weight(grid)
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    major_river_score = np.maximum(micro.river_path_rank, 0.35 * micro.river_rank)
    lake_score = np.maximum(
        micro.lake_basin_rank,
        np.maximum(micro.lake_mask.astype(np.float64), 0.76 * micro.lake_shoreline_rank),
    )
    mouth_score = _river_mouth_zoom_score(grid, refined_rel, micro)
    for label, score in (
        ("major river", major_river_score * lat_weight),
        ("lake basin", lake_score * lat_weight),
        ("river mouth", mouth_score * lat_weight),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=7)
    if not centers:
        centers.append(("hydrology", 0.0, 0.0))
    return centers[:3]


def _fluvial_lacustrine_zoom_centers(
    grid: SphereGrid,
    micro: _Microgeomorphology,
) -> list[tuple[str, float, float]]:
    lat_weight = _zoom_temperate_latitude_weight(grid)
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("meander scrolls", micro.meander_scroll_rank * lat_weight),
        ("lake shoreline", micro.lake_shoreline_rank * lat_weight),
        ("floodplain swales", micro.floodplain_swale_rank * lat_weight),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=14)
    if not centers:
        centers.append(("fluvial / lacustrine microshapes", 0.0, 0.0))
    return centers[:3]


def _zoom_temperate_latitude_weight(
    grid: SphereGrid,
    *,
    preferred_abs_lat: float = 58.0,
    max_abs_lat: float = 72.0,
) -> np.ndarray:
    abs_lat = np.abs(np.asarray(grid.lat, dtype=np.float64))
    span = max(float(max_abs_lat) - float(preferred_abs_lat), 1.0e-9)
    taper = np.clip((float(max_abs_lat) - abs_lat) / span, 0.0, 1.0)
    return np.where(abs_lat < float(max_abs_lat), 0.35 + 0.65 * taper, 0.0)


def _best_zoom_cell(
    grid: SphereGrid,
    score: np.ndarray,
    *,
    blocked: np.ndarray | None,
    min_score: float,
) -> int | None:
    values = np.asarray(score, dtype=np.float64)
    finite = np.isfinite(values)
    eligible = finite & (values > float(min_score))
    if blocked is not None and np.any(eligible & ~blocked):
        eligible &= ~blocked
    if not np.any(eligible):
        return None
    masked = np.where(eligible, values, -np.inf)
    cell = int(np.nanargmax(masked))
    if cell < 0 or cell >= grid.n:
        return None
    return cell


def _river_mouth_zoom_score(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    micro: _Microgeomorphology,
) -> np.ndarray:
    rel = np.asarray(refined_rel, dtype=np.float64)
    land = rel >= 0.0
    coast_distance = np.asarray(micro.land_coast_distance, dtype=np.float64)
    coastal_band = land & (coast_distance >= 0.0) & (coast_distance <= 2.0)
    adjacent_delta = np.zeros(grid.n, dtype=np.float64)
    delta_cells = np.flatnonzero(np.asarray(micro.delta_fan_rank, dtype=np.float64) > 0.0)
    for ocean_cell in delta_cells:
        ocean_cell = int(ocean_cell)
        nbs = grid.neighbors[ocean_cell]
        land_nbs = nbs[land[nbs]]
        if land_nbs.size:
            adjacent_delta[land_nbs] = np.maximum(
                adjacent_delta[land_nbs],
                float(micro.delta_fan_rank[ocean_cell]),
            )
    score = np.zeros(grid.n, dtype=np.float64)
    score[coastal_band] = (
        0.50 * np.asarray(micro.river_path_rank, dtype=np.float64)[coastal_band]
        + 0.22 * np.asarray(micro.river_rank, dtype=np.float64)[coastal_band]
        + 0.32 * np.asarray(micro.delta_plain_mask, dtype=np.float64)[coastal_band]
        + 0.56 * adjacent_delta[coastal_band]
    )
    score[score < 0.04] = 0.0
    return score


def _lowland_alluvial_zoom_centers(
    grid: SphereGrid,
    micro: _Microgeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 62.0
    delta = np.asarray(micro.lowland_alluvial_microrelief_delta, dtype=np.float64)
    centers: list[tuple[str, float, float]] = []
    for label, score in (
        ("alluvial fan", micro.alluvial_fan_rank * low_polar),
        ("lowland plain", micro.lowland_plain_rank * low_polar),
        ("piedmont apron", micro.piedmont_apron_rank * low_polar),
        ("alluvial microrelief", np.abs(delta) * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("lowland / alluvial microrelief", 0.0, 0.0))
    return centers[:3]


def _marine_zoom_centers(
    grid: SphereGrid,
    refined_rel: np.ndarray,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 70.0
    ocean = refined_rel < 0.0
    centers: list[tuple[str, float, float]] = []
    shoal_score = marine.seamount_shoal_mask.astype(np.float64)
    shoal_label = "seamount shoal"
    if np.nanmax(shoal_score) <= 0.0:
        shoal_score = marine.marine_shoal_mask.astype(np.float64)
        shoal_label = "marine shoal"
    for label, score in (
        ("reef morphology", np.maximum.reduce([
            marine.reef_rim_rank,
            marine.atoll_lagoon_rank,
            marine.fringing_reef_rank,
        ]) * low_polar),
        ("islet candidate", marine.island_candidate_rank * low_polar),
        ("reef / atoll", marine.reef_atoll_mask.astype(np.float64) * low_polar),
        (shoal_label, shoal_score * low_polar),
        ("shelf break", marine.shelf_break_rank * low_polar * ocean),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("marine micro", 0.0, 0.0))
    return centers[:3]


def _shelf_slope_microrelief_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    delta = np.asarray(marine.shelf_slope_microrelief_delta, dtype=np.float64)
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("shelf-edge buildup", np.maximum(delta, 0.0) * low_polar),
        ("upper-slope channel", np.maximum(-delta, 0.0) * low_polar),
        ("shelf break", marine.shelf_break_rank * low_polar * (0.35 + np.clip(np.abs(delta), 0.0, 8.0) / 8.0)),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=9)
    if not centers:
        centers.append(("shelf / upper-slope microrelief", 0.0, 0.0))
    return centers[:3]


def _deep_ocean_fabric_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    delta = np.asarray(marine.deep_ocean_fabric_delta, dtype=np.float64)
    centers: list[tuple[str, float, float]] = []
    for label, score in (
        ("fracture-zone trough", marine.fracture_zone_rank * low_polar),
        ("abyssal plain fabric", marine.abyssal_plain_fabric_rank * low_polar),
        ("deep-ocean fabric delta", np.abs(delta) * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("deep-ocean fabric", 0.0, 0.0))
    return centers[:3]


def _reef_atoll_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 62.0
    offshore = _open_ocean_candidate_weight(
        marine.ocean_coast_distance, start=2.0, full=8.0)
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("reef rim", marine.reef_rim_rank * low_polar * (0.52 + 0.48 * offshore)),
        ("atoll lagoon", marine.atoll_lagoon_rank * low_polar * (0.42 + 0.58 * offshore)),
        ("fringing reef", marine.fringing_reef_rank * low_polar),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=6)
    if not centers:
        centers.append(("reef / atoll", 0.0, 0.0))
    return centers[:3]


def _process_island_promotion_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 62.0
    offshore = _open_ocean_candidate_weight(
        marine.ocean_coast_distance, start=2.0, full=8.0)
    promoted_score = (
        marine.process_island_promoted_mask.astype(np.float64)
        * np.maximum(marine.process_island_promotion_rank, 0.75)
        * low_polar
        * (0.56 + 0.44 * offshore)
    )
    atoll_score = (
        marine.atoll_islet_promoted_mask.astype(np.float64)
        * np.maximum(marine.process_island_promotion_rank, 0.75)
        * low_polar
        * (0.46 + 0.54 * offshore)
    )
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("promoted process island", promoted_score),
        ("promoted atoll islet", atoll_score),
        ("islet candidate", marine.island_candidate_rank * low_polar * (0.56 + 0.44 * offshore)),
        ("atoll candidate", marine.atoll_candidate_rank * low_polar * (0.46 + 0.54 * offshore)),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=18)
    if not centers:
        centers.append(("process island promotion", 0.0, 0.0))
    return centers[:3]


def _island_atoll_microshape_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 62.0
    offshore = _open_ocean_candidate_weight(
        marine.ocean_coast_distance, start=2.0, full=8.0)
    promoted_islet = (
        marine.process_island_promoted_mask.astype(np.float64)
        * np.maximum(marine.islet_microshape_rank, 0.75)
        * low_polar
        * (0.56 + 0.44 * offshore)
    )
    promoted_atoll = (
        marine.atoll_islet_promoted_mask.astype(np.float64)
        * np.maximum(marine.atoll_microshape_rank, 0.75)
        * low_polar
        * (0.46 + 0.54 * offshore)
    )
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("subcell promoted islet", promoted_islet),
        ("subcell promoted atoll", promoted_atoll),
        ("subcell islet candidate", marine.islet_microshape_rank * low_polar * (0.56 + 0.44 * offshore)),
        ("subcell atoll candidate", marine.atoll_microshape_rank * low_polar * (0.46 + 0.54 * offshore)),
    ):
        cell = _best_zoom_cell(grid, score, blocked=blocked, min_score=0.0)
        if cell is None:
            continue
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=18)
    if not centers:
        centers.append(("island / atoll microshape", 0.0, 0.0))
    return centers[:3]


def _submarine_highlands_zoom_centers(
    grid: SphereGrid,
    marine: _MarineMicrogeomorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    offshore_gate = np.clip(
        (np.asarray(marine.ocean_coast_distance, dtype=np.float64) - 3.0) / 8.0,
        0.0,
        1.0,
    )
    depth_gate = np.clip(
        (-np.asarray(marine.refined_rel, dtype=np.float64) - 350.0) / 1850.0,
        0.0,
        1.0,
    )
    seamount_score = (
        np.asarray(marine.seamount_peak_rank, dtype=np.float64)
        * low_polar
        * (0.25 + 0.75 * offshore_gate)
        * (0.35 + 0.65 * depth_gate)
    )
    if np.nanmax(seamount_score) <= 0.0:
        seamount_score = marine.seamount_peak_rank * low_polar
    centers: list[tuple[str, float, float]] = []
    for label, score in (
        ("oceanic seamount peak", seamount_score),
        ("oceanic plateau edge", marine.oceanic_plateau_edge_rank * low_polar),
        ("abyssal hill field", marine.abyssal_hill_field_rank * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("submarine highlands", 0.0, 0.0))
    return centers[:3]


def _coastal_zoom_centers(
    grid: SphereGrid,
    coastal: _CoastalMorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    centers: list[tuple[str, float, float]] = []
    for label, score in (
        ("coastal plain", coastal.coastal_plain_rank * low_polar),
        ("coastal cliff", coastal.coastal_cliff_rank * low_polar),
        ("barrier / lagoon", coastal.barrier_lagoon_rank * low_polar),
        ("estuary", coastal.estuary_rank * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("coastal morphology", 0.0, 0.0))
    return centers[:4]


def _coastal_process_linework_zoom_centers(
    grid: SphereGrid,
    coastal: _CoastalMorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    centers: list[tuple[str, float, float]] = []
    blocked = np.zeros(grid.n, dtype=bool)
    for label, score in (
        ("delta distributaries", coastal.delta_distributary_rank * low_polar),
        ("estuary funnel", coastal.estuary_funnel_rank * low_polar),
        ("barrier spit", coastal.barrier_spit_rank * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        selected_score = np.asarray(score, dtype=np.float64).copy()
        selected_score[blocked] = 0.0
        if np.nanmax(selected_score) <= 0.0:
            selected_score = np.asarray(score, dtype=np.float64)
        cell = int(np.nanargmax(selected_score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
        _block_local_neighborhood(grid, blocked, cell, passes=8)
    if not centers:
        centers.append(("coastal process linework", 0.0, 0.0))
    return centers[:3]


def _coastal_depositional_zoom_centers(
    grid: SphereGrid,
    coastal: _CoastalMorphology,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 58.0
    delta = np.asarray(coastal.coastal_depositional_microrelief_delta, dtype=np.float64)
    centers: list[tuple[str, float, float]] = []
    for label, score in (
        ("strandplain", coastal.strandplain_rank * low_polar),
        ("tidal flat", coastal.tidal_flat_rank * low_polar),
        ("depositional plain", coastal.coastal_depositional_plain_rank * low_polar),
        ("coastal depositional microrelief", np.abs(delta) * low_polar),
    ):
        if np.nanmax(score) <= 0.0:
            continue
        cell = int(np.nanargmax(score))
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("coastal depositional microrelief", 0.0, 0.0))
    return centers[:3]


def _render_zoom_sheet(
    grid: SphereGrid,
    parent_rel: np.ndarray,
    refined_rel: np.ndarray,
    detail_delta: np.ndarray,
    path: Path,
    *,
    width: int,
    height: int,
) -> Path:
    parent_r = render.to_raster_continuous(
        grid, parent_rel, width=width, height=height, preserve_sign=True)
    refined_r = render.to_raster_continuous(
        grid, refined_rel, width=width, height=height, preserve_sign=True)
    delta_r = render.to_raster_continuous(
        grid, detail_delta, width=width, height=height)
    centers = _zoom_centers(grid, parent_rel, detail_delta)
    fig, axes = plt.subplots(
        len(centers),
        3,
        figsize=(13.5, 3.9 * len(centers)),
        constrained_layout=True,
    )
    if len(centers) == 1:
        axes = axes[None, :]
    delta_vmax = max(float(np.nanpercentile(np.abs(delta_r), 98.0)), 50.0)
    for row, (label, lon, lat) in enumerate(centers):
        parent_crop, extent = _crop_lon_lat(parent_r, lon, lat)
        refined_crop, _ = _crop_lon_lat(refined_r, lon, lat)
        delta_crop, _ = _crop_lon_lat(delta_r, lon, lat)
        im0 = axes[row, 0].imshow(
            parent_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 0].set_title(f"{label}: parent")
        im1 = axes[row, 1].imshow(
            refined_crop,
            cmap=render.ELEVATION_CMAP,
            norm=render.ELEVATION_NORM,
            extent=extent,
        )
        axes[row, 1].set_title(f"{label}: refined")
        im2 = axes[row, 2].imshow(
            delta_crop,
            cmap="coolwarm",
            vmin=-delta_vmax,
            vmax=delta_vmax,
            extent=extent,
        )
        axes[row, 2].set_title(f"{label}: delta")
        for ax in axes[row, :]:
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
        fig.colorbar(im0, ax=axes[row, 0], shrink=0.75)
        fig.colorbar(im1, ax=axes[row, 1], shrink=0.75)
        fig.colorbar(im2, ax=axes[row, 2], shrink=0.75)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _zoom_centers(
    grid: SphereGrid,
    parent_rel: np.ndarray,
    detail_delta: np.ndarray,
) -> list[tuple[str, float, float]]:
    low_polar = np.abs(grid.lat) < 72.0
    land = (parent_rel >= 0.0) & low_polar
    ocean = (parent_rel < 0.0) & low_polar
    centers: list[tuple[str, float, float]] = []
    for label, mask in (
        ("land detail", land),
        ("ocean detail", ocean),
        ("active relief", low_polar),
    ):
        if not mask.any():
            continue
        values = np.abs(detail_delta)
        cell = int(np.flatnonzero(mask)[np.argmax(values[mask])])
        centers.append((label, float(grid.lon[cell]), float(grid.lat[cell])))
    if not centers:
        centers.append(("global detail", 0.0, 0.0))
    return centers[:3]


def _crop_lon_lat(
    raster: np.ndarray,
    lon: float,
    lat: float,
    *,
    lon_width: float = 72.0,
    lat_height: float = 36.0,
) -> tuple[np.ndarray, list[float]]:
    h, w = raster.shape
    lon_axis = np.linspace(-180.0, 180.0, w)
    lat_axis = np.linspace(90.0, -90.0, h)
    lon_min = max(-180.0, float(lon) - lon_width * 0.5)
    lon_max = min(180.0, float(lon) + lon_width * 0.5)
    lat_min = max(-90.0, float(lat) - lat_height * 0.5)
    lat_max = min(90.0, float(lat) + lat_height * 0.5)
    cols = np.where((lon_axis >= lon_min) & (lon_axis <= lon_max))[0]
    rows = np.where((lat_axis >= lat_min) & (lat_axis <= lat_max))[0]
    if cols.size == 0:
        cols = np.asarray([int(np.argmin(np.abs(lon_axis - lon)))])
    if rows.size == 0:
        rows = np.asarray([int(np.argmin(np.abs(lat_axis - lat)))])
    crop = raster[rows.min():rows.max() + 1, cols.min():cols.max() + 1]
    extent = [
        float(lon_axis[cols.min()]),
        float(lon_axis[cols.max()]),
        float(lat_axis[rows.max()]),
        float(lat_axis[rows.min()]),
    ]
    return crop, extent


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)
