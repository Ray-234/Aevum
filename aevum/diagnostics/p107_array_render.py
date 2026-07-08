"""Render P107 terminal audit arrays without re-running a world.

The P107 audit stores enough terminal fields in ``p107_terminal_arrays.npz`` to
inspect the high-resolution map after a no-render run.  This helper renders a
fixed visual QA pack from those arrays so expensive 24000-cell simulations do
not need to be repeated just to produce PNGs.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

from aevum import render
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS


SCHEMA = "aevum.p107_array_render.v1"


def render_p107_array_assets(
    source: str | Path,
    outdir: str | Path,
    *,
    width: int = 720,
    height: int = 360,
) -> dict[str, Any]:
    """Render PNG assets from a P107 run directory or terminal metrics file."""
    metrics_path = _resolve_metrics_path(Path(source))
    metrics = json.loads(metrics_path.read_text())
    arrays_path = _resolve_arrays_path(metrics_path, metrics)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    with np.load(arrays_path) as archive:
        arrays = {key: archive[key] for key in archive.files}

    n = int(arrays["grid_lat"].shape[0])
    grid = SphereGrid.fibonacci(n, CONSTANTS.EARTH_RADIUS)
    idx = _raster_index(grid, width=width, height=height)
    sea_level = float(metrics.get("sea_level_m", 0.0))
    manifest = metrics.get("array_archive", {}).get("manifest", {})
    fields = manifest.get("fields", {})

    assets: dict[str, str] = {}

    elev = _array_field(arrays, fields, "terrain.elevation_m").astype(np.float64)
    rel = elev - sea_level
    assets["elevation.png"] = str(_render_elevation(rel, grid, idx, out))
    assets["bathymetry_shelf_slope_abyss.png"] = str(
        _render_bathymetry(rel, grid, idx, out))

    plate = _array_field(arrays, fields, "tectonics.plate_id").astype(np.float64)
    assets["plates.png"] = str(_render_discrete(
        plate, grid, idx, out / "plates.png", "Plates", cmap="tab20"))

    crust_type = _array_field(arrays, fields, "crust.type").astype(np.float64)
    crust_age = _array_field(arrays, fields, "crust.age_myr").astype(np.float64)
    assets["crust_age.png"] = str(
        _render_crust_age(crust_age, crust_type, grid, idx, out))

    ocean_basin = _array_field(arrays, fields, "ocean.basin_id", required=False)
    if ocean_basin is not None:
        assets["ocean_basins.png"] = str(_render_discrete(
            ocean_basin.astype(np.float64),
            grid,
            idx,
            out / "ocean_basins.png",
            "Ocean basins",
            cmap="tab20",
            vmin=-1.0,
        ))

    depth_province = _array_field(
        arrays, fields, "ocean.depth_province", required=False)
    if depth_province is not None:
        assets["ocean_depth_provinces.png"] = str(_render_ocean_depth_provinces(
            depth_province.astype(np.float64), grid, idx, out))

    terrain_province = _array_field(
        arrays, fields, "terrain.province", required=False)
    if terrain_province is not None:
        assets["terrain_provinces.png"] = str(_render_discrete(
            terrain_province.astype(np.float64),
            grid,
            idx,
            out / "terrain_provinces.png",
            "Terrain provinces",
            cmap="tab20",
            vmin=0.0,
        ))

    detail = _array_field(
        arrays, fields, "terrain.continental_detail", required=False)
    detail_region = _array_field(
        arrays, fields, "terrain.continental_detail_region_code",
        required=False)
    if detail is not None:
        detail = detail.astype(np.float64)
        assets["continental_detail_raw_provinces.png"] = str(
            _render_continental_detail_codes(
                detail,
                grid,
                idx,
                out / "continental_detail_raw_provinces.png",
                "P107 array continental detail raw provinces",
            )
        )
        if detail_region is None:
            assets["continental_detail_provinces.png"] = str(
                _render_continental_detail_codes(
                    detail,
                    grid,
                    idx,
                    out / "continental_detail_provinces.png",
                    "P107 array continental detail provinces",
                )
            )
    if detail_region is not None:
        detail_region = detail_region.astype(np.float64)
        assets["continental_detail_region_provinces.png"] = str(
            _render_continental_detail_codes(
                detail_region,
                grid,
                idx,
                out / "continental_detail_region_provinces.png",
                "P107 array continental detail regional expression",
            )
        )
        assets["continental_detail_provinces.png"] = str(
            _render_continental_detail_codes(
                detail_region,
                grid,
                idx,
                out / "continental_detail_provinces.png",
                "P107 array continental detail provinces",
            )
        )

    internal_block_region = _array_field(
        arrays, fields, "terrain.internal_geographic_block_region_code",
        required=False)
    if internal_block_region is not None:
        assets["internal_geographic_block_region_code.png"] = str(
            _render_internal_block_codes(
                internal_block_region.astype(np.float64),
                grid,
                idx,
                out / "internal_geographic_block_region_code.png",
            )
        )

    inland_region = _array_field(
        arrays, fields, "terrain.inland_geomorphology_region_code",
        required=False)
    if inland_region is not None:
        assets["inland_geomorphology_regions.png"] = str(
            _render_inland_geomorphology_codes(
                inland_region.astype(np.float64),
                grid,
                idx,
                out / "inland_geomorphology_regions.png",
            )
        )

    terrain_cont_province = _array_field(
        arrays, fields, "terrain.continental_province_code", required=False)
    if terrain_cont_province is not None:
        assets["terrain_continental_province_code.png"] = str(_render_discrete(
            terrain_cont_province.astype(np.float64),
            grid,
            idx,
            out / "terrain_continental_province_code.png",
            "Terrain continental province code",
            cmap="tab20",
            vmin=0.0,
        ))

    orogenic_hierarchy = _array_field(
        arrays, fields, "terrain.orogenic_parent_hierarchy", required=False)
    orogenic_hierarchy_raster = None
    if orogenic_hierarchy is not None:
        orogenic_hierarchy = orogenic_hierarchy.astype(np.float64)
        orogenic_hierarchy_raster = _raster(orogenic_hierarchy, idx)
        assets["orogenic_parent_hierarchy.png"] = str(
            _render_orogenic_parent_hierarchy(
                orogenic_hierarchy, grid, idx, out))
        assets["orogenic_parent_hierarchy_overlay.png"] = str(
            _render_orogenic_parent_hierarchy_overlay(
                orogenic_hierarchy, rel, grid, idx, out))

    orogenic_spine = _array_field(
        arrays, fields, "terrain.orogenic_hierarchy_spine", required=False)
    orogenic_spine_raster = None
    if orogenic_spine is not None:
        orogenic_spine = orogenic_spine.astype(np.float64)
        orogenic_spine_raster = _raster(orogenic_spine, idx)
        assets["orogenic_hierarchy_spines.png"] = str(
            _render_orogenic_hierarchy_spines(orogenic_spine, grid, idx, out))
        assets["orogenic_hierarchy_spine_overlay.png"] = str(
            _render_orogenic_hierarchy_spine_overlay(
                orogenic_spine, rel, grid, idx, out))

    orogenic_halo = _array_field(
        arrays, fields, "terrain.orogenic_shoulder_halo", required=False)
    orogenic_halo_raster = None
    if orogenic_halo is not None:
        orogenic_halo = orogenic_halo.astype(np.float64)
        orogenic_halo_raster = _raster(orogenic_halo, idx)
        assets["orogenic_shoulder_halo.png"] = str(
            _render_orogenic_shoulder_halo(orogenic_halo, grid, idx, out))

    orogenic_apron = _array_field(
        arrays, fields, "terrain.orogenic_highland_apron", required=False)
    orogenic_apron_raster = None
    if orogenic_apron is not None:
        orogenic_apron = orogenic_apron.astype(np.float64)
        orogenic_apron_raster = _raster(orogenic_apron, idx)
        assets["orogenic_highland_apron.png"] = str(
            _render_orogenic_highland_apron(orogenic_apron, grid, idx, out))

    if (
        orogenic_hierarchy is not None
        or orogenic_spine is not None
        or orogenic_halo is not None
        or orogenic_apron is not None
    ):
        assets["orogenic_belt_morphology_overlay.png"] = str(
            _render_orogenic_belt_morphology_overlay(
                orogenic_hierarchy,
                orogenic_spine,
                orogenic_halo,
                orogenic_apron,
                rel,
                grid,
                idx,
                out,
            )
        )

    boundary = _combined_mask_raster(
        arrays,
        idx,
        [
            ("boundary__ridge", 1),
            ("boundary__transform", 2),
            ("boundary__trench", 3),
            ("boundary__suture", 4),
            ("boundary__active_margin", 5),
            ("boundary__passive_margin", 6),
            ("boundary__collision", 7),
        ],
    )
    assets["tectonic_boundaries.png"] = str(
        _render_boundary_raster(boundary, out / "tectonic_boundaries.png"))

    object_raster = _combined_mask_raster(
        arrays,
        idx,
        [
            ("object__province_mid_ocean_ridge", 1),
            ("object__fracture_zone", 2),
            ("object__margin_trench", 3),
            ("object__island_arc", 4),
            ("object__back_arc_basin", 5),
            ("object__seamount_chain", 6),
            ("object__oceanic_plateau", 7),
            ("object__microcontinent", 8),
            ("object__mountain_orogen", 9),
            ("object__mountain_old_orogen", 10),
            ("object__mountain_plateau", 11),
            ("object__parent_orogen_highland_apron", 12),
        ],
    )
    assets["object_masks.png"] = str(
        _render_object_raster(object_raster, out / "object_masks.png"))

    contact_panels = {
        "Elevation": _raster_continuous(rel, grid, idx, preserve_sign=True),
        "Bathymetry": _bathymetry_raster(rel, grid, idx),
        "Plates": _raster(plate, idx),
        "Crust age": _raster(crust_age, idx),
    }
    if orogenic_hierarchy_raster is not None:
        contact_panels["Orogen hierarchy"] = orogenic_hierarchy_raster
    if orogenic_spine_raster is not None:
        contact_panels["Orogen spines"] = orogenic_spine_raster
    if orogenic_halo_raster is not None:
        contact_panels["Orogen halo"] = orogenic_halo_raster
    if orogenic_apron_raster is not None:
        contact_panels["Orogen apron"] = orogenic_apron_raster
    contact_panels["Boundaries"] = boundary
    contact_panels["Objects"] = object_raster
    assets["p107_array_contact_sheet.png"] = str(
        _render_contact_sheet(
            contact_panels,
            out / "p107_array_contact_sheet.png",
        )
    )

    summary = {
        "schema": SCHEMA,
        "source_metrics": str(metrics_path),
        "source_arrays": str(arrays_path),
        "cells": int(n),
        "width": int(width),
        "height": int(height),
        "sea_level_m": sea_level,
        "assets": assets,
        "p110a_summary": metrics.get("p110a_modern_planform", {}).get(
            "summary", {}),
        "p110a_warning_flags": metrics.get(
            "p110a_modern_planform", {}).get("warning_flags", []),
        "p109_pass": metrics.get(
            "p109_hypsometry_comparison", {}).get("within_p109_envelope"),
    }
    summary_path = out / "p107_array_render_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=_json_default) + "\n")
    return summary


def _resolve_metrics_path(source: Path) -> Path:
    source = source.expanduser()
    if source.is_dir():
        candidate = source / "p107_terminal_metrics.json"
    else:
        candidate = source
    if not candidate.exists():
        raise FileNotFoundError(f"P107 terminal metrics not found: {candidate}")
    return candidate


def _resolve_arrays_path(metrics_path: Path, metrics: dict[str, Any]) -> Path:
    raw = str(metrics.get("array_archive", {}).get("path", ""))
    if not raw:
        raise ValueError("metrics JSON does not contain array_archive.path")
    path = Path(raw).expanduser()
    if path.exists():
        return path
    joined = metrics_path.parent / path
    if joined.exists():
        return joined
    sibling = metrics_path.parent / "p107_terminal_arrays.npz"
    if sibling.exists():
        return sibling
    raise FileNotFoundError(f"P107 terminal arrays not found: {raw}")


def _raster_index(grid: SphereGrid, *, width: int, height: int) -> np.ndarray:
    lon = np.linspace(-180.0, 180.0, int(width))
    lat = np.linspace(90.0, -90.0, int(height))
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")
    return grid.nearest_latlon(lat_grid.ravel(), lon_grid.ravel()).reshape(
        int(height), int(width))


def _raster(values: np.ndarray, idx: np.ndarray) -> np.ndarray:
    return np.asarray(values)[idx]


def _raster_continuous(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    *,
    preserve_sign: bool = False,
) -> np.ndarray:
    return render.to_raster_continuous(
        grid,
        np.asarray(values),
        width=int(idx.shape[1]),
        height=int(idx.shape[0]),
        preserve_sign=preserve_sign,
    )


def _array_field(
    arrays: dict[str, np.ndarray],
    fields: dict[str, str],
    field_name: str,
    *,
    required: bool = True,
) -> np.ndarray | None:
    key = fields.get(field_name)
    if key is None:
        key = "field__" + field_name.replace(".", "_")
    if key not in arrays:
        if required:
            raise KeyError(f"missing P107 array field {field_name!r} ({key})")
        return None
    return arrays[key]


def _render_elevation(
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    raster = _raster_continuous(rel, grid, idx, preserve_sign=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = render.render_elevation_raster(
        ax,
        raster,
        title=f"P107 array elevation (m rel. sea level), {grid.n} cells",
    )
    render.add_elevation_colorbar(fig, ax, im)
    path = outdir / "elevation.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _bathymetry_raster(
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
) -> np.ndarray:
    depth = np.where(np.asarray(rel) < 0.0, -np.asarray(rel), np.nan)
    return _raster_continuous(depth, grid, idx)


def _render_bathymetry(
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    raster = np.ma.masked_invalid(_bathymetry_raster(rel, grid, idx))
    values = -np.asarray(rel)[np.asarray(rel) < 0.0]
    vmax = max(float(np.nanpercentile(values, 98.0)) if values.size else 1.0, 1000.0)
    cmap = plt.get_cmap("cividis_r").with_extremes(bad="#efe8d7")
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        raster,
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(f"P107 array bathymetry: shelf, slope, abyss, {grid.n} cells")
    fig.colorbar(im, ax=ax, shrink=0.7)
    path = outdir / "bathymetry_shelf_slope_abyss.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_discrete(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    path: Path,
    title: str,
    *,
    cmap: str = "tab20",
    vmin: float | None = None,
) -> Path:
    raster = _raster(values, idx)
    vmax = max(float(np.nanmax(values)), 1.0)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(
        raster,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(f"P107 array {title.lower()}, {grid.n} cells")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_crust_age(
    age: np.ndarray,
    ctype: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    age_r = _raster(age, idx)
    type_r = _raster(ctype, idx)
    ocean = type_r < 0.5
    continent = ~ocean
    ocean_age = np.ma.masked_where(~ocean, age_r)
    cont_age = np.ma.masked_where(~continent, age_r)
    ocean_values = age[ctype < 0.5]
    cont_values = age[ctype >= 0.5]
    ocean_vmax = 300.0
    if ocean_values.size:
        ocean_vmax = float(np.clip(np.percentile(ocean_values, 98.0), 120.0, 500.0))
    cont_vmax = 3000.0
    if cont_values.size:
        cont_vmax = float(np.clip(np.percentile(cont_values, 98.0), 900.0, 3500.0))
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_facecolor("#efe8d7")
    ax.imshow(
        cont_age,
        cmap=render.CONTINENT_CRUST_AGE_CMAP,
        vmin=0.0,
        vmax=cont_vmax,
        extent=[-180, 180, -90, 90],
        alpha=0.58,
    )
    im = ax.imshow(
        ocean_age,
        cmap=render.OCEAN_CRUST_AGE_CMAP,
        vmin=0.0,
        vmax=ocean_vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(f"P107 array oceanic crust age (Myr), {grid.n} cells")
    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label("oceanic crust age (Myr)")
    path = outdir / "crust_age.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_ocean_depth_provinces(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    colors = [
        "#efe8d7",
        "#9bd4d8",
        "#4aa6b5",
        "#287c9d",
        "#173b57",
        "#c7b86a",
        "#301934",
        "#d98ca4",
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    raster = _raster(values, idx)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title(f"P107 array ocean depth provinces, {grid.n} cells")
    fig.colorbar(im, ax=ax, shrink=0.7)
    path = outdir / "ocean_depth_provinces.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_continental_detail_codes(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    path: Path,
    title: str,
) -> Path:
    colors = [
        "#173b57",
        "#6f8e59",
        "#4f9a52",
        "#b7c88a",
        "#87b6a6",
        "#b35d4d",
        "#8f6fb2",
        "#c99454",
    ]
    labels = [
        "ocean",
        "shield",
        "platform",
        "basin",
        "rift/passive basin",
        "orogen",
        "plateau",
        "arc microcontinent",
    ]
    return _render_coded_raster(_raster(values, idx), path, title, colors, labels)


def _render_internal_block_codes(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    path: Path,
) -> Path:
    colors = [
        "#173b57",
        "#6f8e59",
        "#4f9a52",
        "#b7c88a",
        "#87b6a6",
        "#b35d4d",
        "#c99454",
    ]
    labels = [
        "none/ocean",
        "craton core",
        "stable platform",
        "intracratonic basin",
        "mobile belt",
        "rifted margin",
        "accreted terrane",
    ]
    return _render_coded_raster(
        _raster(values, idx),
        path,
        f"P107 array internal geographic block regional expression, {grid.n} cells",
        colors,
        labels,
    )


def _render_inland_geomorphology_codes(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    path: Path,
) -> Path:
    colors = [
        "#173b57",
        "#6f8e59",
        "#5f9f4a",
        "#b7c88a",
        "#b35d4d",
        "#87b6a6",
        "#c7aa58",
        "#d98a34",
        "#8f6fb2",
    ]
    labels = [
        "none/ocean",
        "shield",
        "platform",
        "sag basin",
        "old orogen",
        "rift",
        "platform swell",
        "escarpment",
        "plateau margin",
    ]
    return _render_coded_raster(
        _raster(values, idx),
        path,
        f"P107 array inland geomorphology regions, {grid.n} cells",
        colors,
        labels,
    )


def _combined_mask_raster(
    arrays: dict[str, np.ndarray],
    idx: np.ndarray,
    specs: list[tuple[str, int]],
) -> np.ndarray:
    field = np.zeros(next(iter(arrays.values())).shape[0], dtype=np.float64)
    for key, code in specs:
        mask = arrays.get(key)
        if mask is None:
            continue
        field[np.asarray(mask, dtype=bool)] = float(code)
    return _raster(field, idx)


def _render_boundary_raster(raster: np.ndarray, path: Path) -> Path:
    colors = [
        "#f7f7f7",
        "#d73027",
        "#4575b4",
        "#542788",
        "#7b3294",
        "#f46d43",
        "#66bd63",
        "#000000",
    ]
    labels = [
        "none",
        "ridge",
        "transform",
        "trench",
        "suture",
        "active margin",
        "passive margin",
        "collision",
    ]
    return _render_coded_raster(raster, path, "P107 array tectonic boundaries",
                                colors, labels)


def _render_object_raster(raster: np.ndarray, path: Path) -> Path:
    colors = [
        "#f7f7f7",
        "#d73027",
        "#4575b4",
        "#542788",
        "#fdae61",
        "#abd9e9",
        "#1a9850",
        "#fee08b",
        "#984ea3",
        "#000000",
        "#999999",
        "#e7298a",
        "#c2a5cf",
    ]
    labels = [
        "none",
        "ridge province",
        "fracture zone",
        "trench",
        "island arc",
        "back-arc",
        "seamount chain",
        "oceanic plateau",
        "microcontinent",
        "orogen",
        "old orogen",
        "plateau",
        "orogen apron",
    ]
    return _render_coded_raster(raster, path, "P107 array tectonic object masks",
                                colors, labels)


def _render_orogenic_parent_hierarchy(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    colors = [
        "#f7f7f7",
        "#b8d88e",
        "#d98a34",
        "#6f3f8f",
    ]
    labels = [
        "none",
        "foreland slope",
        "branch range",
        "main crest",
    ]
    raster = _raster(values, idx)
    return _render_coded_raster(
        raster,
        outdir / "orogenic_parent_hierarchy.png",
        f"P107 array orogenic parent hierarchy, {grid.n} cells",
        colors,
        labels,
    )


def _render_orogenic_parent_hierarchy_overlay(
    values: np.ndarray,
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    hierarchy_raster = _raster(values, idx)
    elevation_raster = _raster(rel, idx)
    masked = np.ma.masked_where(hierarchy_raster < 0.5, hierarchy_raster)
    cmap = ListedColormap(["#b8d88e", "#d98a34", "#6f3f8f"])
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], cmap.N)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elevation_raster,
        title=f"P107 array orogenic hierarchy on elevation, {grid.n} cells",
    )
    im = ax.imshow(
        masked,
        cmap=cmap,
        norm=norm,
        alpha=0.74,
        extent=[-180, 180, -90, 90],
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=[1, 2, 3])
    cb.ax.set_yticklabels(["foreland slope", "branch range", "main crest"])
    path = outdir / "orogenic_parent_hierarchy_overlay.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_orogenic_hierarchy_spines(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    colors = [
        "#f7f7f7",
        "#d98a34",
        "#6f3f8f",
    ]
    labels = [
        "none",
        "branch spine",
        "main crest spine",
    ]
    raster = _raster(np.where(values >= 3.0, 2.0, np.where(values >= 2.0, 1.0, 0.0)), idx)
    return _render_coded_raster(
        raster,
        outdir / "orogenic_hierarchy_spines.png",
        f"P107 array orogenic hierarchy spines, {grid.n} cells",
        colors,
        labels,
    )


def _render_orogenic_hierarchy_spine_overlay(
    values: np.ndarray,
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    spine_raster = _raster(
        np.where(values >= 3.0, 2.0, np.where(values >= 2.0, 1.0, 0.0)),
        idx,
    )
    elevation_raster = _raster(rel, idx)
    masked = np.ma.masked_where(spine_raster < 0.5, spine_raster)
    cmap = ListedColormap(["#d98a34", "#6f3f8f"])
    norm = BoundaryNorm([0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elevation_raster,
        title=f"P107 array orogenic spines on elevation, {grid.n} cells",
    )
    im = ax.imshow(
        masked,
        cmap=cmap,
        norm=norm,
        alpha=0.92,
        extent=[-180, 180, -90, 90],
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=[1, 2])
    cb.ax.set_yticklabels(["branch spine", "main crest spine"])
    path = outdir / "orogenic_hierarchy_spine_overlay.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_orogenic_shoulder_halo(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    colors = [
        "#f7f7f7",
        "#66c2a5",
    ]
    labels = [
        "none",
        "shoulder/foothill halo",
    ]
    raster = _raster(np.where(values > 0.0, 1.0, 0.0), idx)
    return _render_coded_raster(
        raster,
        outdir / "orogenic_shoulder_halo.png",
        f"P107 array orogenic shoulder halo, {grid.n} cells",
        colors,
        labels,
    )


def _render_orogenic_highland_apron(
    values: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    colors = [
        "#f7f7f7",
        "#c2a5cf",
    ]
    labels = [
        "none",
        "highland apron",
    ]
    raster = _raster(np.where(values > 0.0, 1.0, 0.0), idx)
    return _render_coded_raster(
        raster,
        outdir / "orogenic_highland_apron.png",
        f"P107 array orogenic highland apron, {grid.n} cells",
        colors,
        labels,
    )


def _render_orogenic_belt_morphology_overlay(
    hierarchy: np.ndarray | None,
    spine: np.ndarray | None,
    halo: np.ndarray | None,
    apron: np.ndarray | None,
    rel: np.ndarray,
    grid: SphereGrid,
    idx: np.ndarray,
    outdir: Path,
) -> Path:
    code = np.zeros(grid.n, dtype=np.float64)
    if hierarchy is not None:
        hierarchy = np.asarray(hierarchy, dtype=np.float64)
        code[hierarchy >= 1.0] = 3.0
        code[hierarchy >= 2.0] = 4.0
        code[hierarchy >= 3.0] = 5.0
    if halo is not None:
        code[np.asarray(halo, dtype=np.float64) > 0.0] = 1.0
    if apron is not None:
        code[np.asarray(apron, dtype=np.float64) > 0.0] = 2.0
    if spine is not None:
        spine = np.asarray(spine, dtype=np.float64)
        code[spine >= 2.0] = 6.0
        code[spine >= 3.0] = 7.0
    code_raster = _raster(code, idx)
    elevation_raster = _raster(rel, idx)
    masked = np.ma.masked_where(code_raster < 0.5, code_raster)
    cmap = ListedColormap([
        "#66c2a5",
        "#c2a5cf",
        "#b8d88e",
        "#d98a34",
        "#6f3f8f",
        "#ff8f00",
        "#d50000",
    ])
    norm = BoundaryNorm(
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5],
        cmap.N,
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    render.render_elevation_raster(
        ax,
        elevation_raster,
        title=f"P107 array orogenic belt morphology on elevation, {grid.n} cells",
    )
    im = ax.imshow(
        masked,
        cmap=cmap,
        norm=norm,
        alpha=0.76,
        extent=[-180, 180, -90, 90],
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=[1, 2, 3, 4, 5, 6, 7])
    cb.ax.set_yticklabels([
        "shoulder halo",
        "highland apron",
        "foreland slope",
        "branch range",
        "main crest",
        "branch spine",
        "crest spine",
    ])
    path = outdir / "orogenic_belt_morphology_overlay.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_coded_raster(
    raster: np.ndarray,
    path: Path,
    title: str,
    colors: list[str],
    labels: list[str],
) -> Path:
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=np.arange(len(labels)))
    cb.ax.set_yticklabels(labels)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_contact_sheet(rasters: dict[str, np.ndarray], path: Path) -> Path:
    cols = 2
    rows = max(1, (len(rasters) + cols - 1) // cols)
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(14, max(4.5, 3.0 * rows)),
        constrained_layout=True,
    )
    axes_flat = np.asarray(axes).reshape(-1)
    for ax, (title, raster) in zip(axes_flat, rasters.items()):
        if title == "Elevation":
            im = render.render_elevation_raster(ax, raster, title=title)
            fig.colorbar(im, ax=ax, shrink=0.68)
        elif title == "Bathymetry":
            masked = np.ma.masked_invalid(raster)
            im = ax.imshow(masked, cmap="cividis_r", extent=[-180, 180, -90, 90])
            ax.set_title(title)
            fig.colorbar(im, ax=ax, shrink=0.68)
        elif title == "Crust age":
            im = ax.imshow(raster, cmap=render.OCEAN_CRUST_AGE_CMAP,
                           extent=[-180, 180, -90, 90])
            ax.set_title(title)
            fig.colorbar(im, ax=ax, shrink=0.68)
        elif title == "Boundaries":
            cmap = ListedColormap([
                "#f7f7f7", "#d73027", "#4575b4", "#542788",
                "#7b3294", "#f46d43", "#66bd63", "#000000",
            ])
            norm = BoundaryNorm(np.arange(-0.5, 8.5), cmap.N)
            ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
            ax.set_title(title)
        elif title == "Objects":
            cmap = ListedColormap([
                "#f7f7f7", "#d73027", "#4575b4", "#542788",
                "#fdae61", "#abd9e9", "#1a9850", "#fee08b",
                "#984ea3", "#000000", "#999999", "#e7298a",
                "#c2a5cf",
            ])
            norm = BoundaryNorm(np.arange(-0.5, 13.5), cmap.N)
            ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
            ax.set_title(title)
        elif title == "Orogen hierarchy":
            cmap = ListedColormap([
                "#f7f7f7", "#b8d88e", "#d98a34", "#6f3f8f",
            ])
            norm = BoundaryNorm(np.arange(-0.5, 4.5), cmap.N)
            ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
            ax.set_title(title)
        elif title == "Orogen spines":
            cmap = ListedColormap([
                "#f7f7f7", "#d98a34", "#6f3f8f",
            ])
            norm = BoundaryNorm(np.arange(-0.5, 3.5), cmap.N)
            ax.imshow(
                np.where(raster >= 3.0, 2.0, np.where(raster >= 2.0, 1.0, 0.0)),
                cmap=cmap,
                norm=norm,
                extent=[-180, 180, -90, 90],
            )
            ax.set_title(title)
        elif title == "Orogen halo":
            cmap = ListedColormap(["#f7f7f7", "#66c2a5"])
            norm = BoundaryNorm(np.arange(-0.5, 2.5), cmap.N)
            ax.imshow(
                np.where(raster > 0.0, 1.0, 0.0),
                cmap=cmap,
                norm=norm,
                extent=[-180, 180, -90, 90],
            )
            ax.set_title(title)
        elif title == "Orogen apron":
            cmap = ListedColormap(["#f7f7f7", "#c2a5cf"])
            norm = BoundaryNorm(np.arange(-0.5, 2.5), cmap.N)
            ax.imshow(
                np.where(raster > 0.0, 1.0, 0.0),
                cmap=cmap,
                norm=norm,
                extent=[-180, 180, -90, 90],
            )
            ax.set_title(title)
        else:
            im = ax.imshow(raster, cmap="tab20", extent=[-180, 180, -90, 90])
            ax.set_title(title)
            fig.colorbar(im, ax=ax, shrink=0.68)
    for ax in axes_flat[len(rasters):]:
        ax.axis("off")
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)
