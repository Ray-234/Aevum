"""Real-Earth distribution screening for tectonic geomorphology.

This module is an ETOPO-anchored envelope check, not a full ETOPO raster fit.
It turns the current world's hypsometry, ocean provinces, morphology, and
object coverage into comparable scalar diagnostics so later tuning can move
from visual inspection to repeatable distribution checks.
"""
from __future__ import annotations

from typing import Any

import numpy as np


REFERENCE_SOURCES = [
    {
        "id": "NOAA_ETOPO_2022",
        "title": "NOAA NCEI ETOPO 2022 Global Relief Model",
        "url": "https://www.ncei.noaa.gov/products/etopo-global-relief-model",
        "doi": "10.25921/fd45-gt74",
        "role": "global topography and bathymetry reference for later raster-derived envelopes",
    },
    {
        "id": "Aevum_ETOP0_envelope_stage",
        "title": "Aevum initial real-Earth distribution envelopes",
        "url": "docs/EARTH_GEOMORPHOLOGY_COVERAGE.md",
        "role": "screening ranges pending direct ETOPO raster sampling",
    },
]


REFERENCE_ENVELOPES: dict[str, dict[str, Any]] = {
    "land_fraction": {
        "range": (0.25, 0.33),
        "basis": "Modern Earth has about 29 percent exposed land.",
    },
    "largest_land_component_fraction": {
        "range": (0.20, 0.60),
        "basis": "Modern Earth has several large continents, not one dominant supercontinent.",
    },
    "land_component_count": {
        "range": (4.0, 14.0),
        "basis": "Coarse Earth-like review should retain multiple major land components.",
    },
    "land_ribbon_fraction_gt_0_5": {
        "range": (0.0, 0.35),
        "basis": "Mainland area should not be dominated by one-cell-wide or chain-like land.",
    },
    "land_coastline_complexity_largest": {
        "range": (0.0, 8.0),
        "basis": "Largest landmass should not have low-resolution coastline over-complexity.",
    },
    "land_elevation_mean_m": {
        "range": (400.0, 1200.0),
        "basis": "Modern continental mean elevation is roughly order 0.8 km above sea level.",
    },
    "land_elevation_p95_m": {
        "range": (1500.0, 5200.0),
        "basis": "Earth-like worlds need a high-relief tail but not all land as highland.",
    },
    "high_land_fraction_gt2500m": {
        "range": (0.0, 0.18),
        "basis": "High plateaus and mountain belts should be limited-area provinces.",
    },
    "orogen_or_plateau_fraction_of_land": {
        "range": (0.02, 0.35),
        "basis": "Orogens and plateaus should be present but not overpaint continents.",
    },
    "shelf_fraction_of_ocean": {
        "range": (0.04, 0.14),
        "basis": "Continental shelves are a small but visible fraction of global ocean.",
    },
    "slope_rise_fraction_of_ocean": {
        "range": (0.10, 0.35),
        "basis": "Continental slopes/rises should bridge shelves and abyssal plains.",
    },
    "abyss_fraction_of_ocean": {
        "range": (0.30, 0.65),
        "basis": "A mature ocean should have broad abyssal-plain/deep-basin area.",
    },
    "ridge_fraction_of_ocean": {
        "range": (0.02, 0.25),
        "basis": "Ridges should be globally connected but not dominate the ocean floor.",
    },
    "trench_fraction_of_ocean": {
        "range": (0.005, 0.12),
        "basis": "Trenches are narrow active-margin features, not broad provinces.",
    },
    "shelf_to_abyss_depth_delta_m": {
        "range": (1200.0, 5200.0),
        "basis": "Shelf seas and abyssal plains must be bathymetrically separated.",
    },
    "nearshore_superdeep_fraction_gt2500m": {
        "range": (0.0, 0.03),
        "basis": "Superdeep water immediately offshore should be rare outside trenches.",
    },
    "far_ocean_shallow_fraction_lt1500m": {
        "range": (0.0, 0.12),
        "basis": "Remote open ocean should mostly be deep basin, not shallow shelf.",
    },
}


def earth_reference_distribution_metrics(
    world,
    tectonics: dict[str, Any] | None = None,
    geomorphology: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare a world with initial real-Earth distribution envelopes."""
    if hasattr(world, "world"):
        engine = world
        world = engine.world
        if tectonics is None:
            from aevum import validation

            tectonics = validation.tectonic_diagnostics(engine)
    if tectonics is None:
        raise ValueError("tectonics diagnostics must be supplied when passing a WorldState")

    grid = world.grid
    area = grid.cell_area
    total_area = float(area.sum())
    rel = np.asarray(world.get_field("terrain.elevation_m", 0.0), dtype=np.float64) - float(world.sea_level)
    land = rel >= 0.0

    morphology = tectonics.get("morphology", {})
    exposed = morphology.get("exposed_land", {})
    ocean = tectonics.get("ocean_geography", {})

    values: dict[str, float] = {
        "land_fraction": float(area[land].sum() / total_area) if total_area else 0.0,
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
        "shelf_fraction_of_ocean": float(ocean.get("shelf_fraction_of_ocean", 0.0)),
        "slope_rise_fraction_of_ocean": float(ocean.get("slope_rise_fraction_of_ocean", 0.0)),
        "abyss_fraction_of_ocean": float(ocean.get("abyss_fraction_of_ocean", 0.0)),
        "ridge_fraction_of_ocean": float(ocean.get("ridge_fraction_of_ocean", 0.0)),
        "trench_fraction_of_ocean": float(ocean.get("trench_fraction_of_ocean", 0.0)),
        "shelf_to_abyss_depth_delta_m": float(ocean.get("shelf_to_abyss_depth_delta_m", 0.0)),
        "nearshore_superdeep_fraction_gt2500m": float(
            ocean.get("nearshore_superdeep_fraction_gt2500m", 0.0)
        ),
        "far_ocean_shallow_fraction_lt1500m": float(
            ocean.get("far_ocean_shallow_fraction_lt1500m", 0.0)
        ),
    }

    if land.any():
        land_rel = rel[land]
        values["land_elevation_mean_m"] = float(np.average(land_rel, weights=area[land]))
        values["land_elevation_p95_m"] = float(np.percentile(land_rel, 95))
        values["high_land_fraction_gt2500m"] = float(
            area[land & (rel > 2500.0)].sum() / max(area[land].sum(), 1.0)
        )
    else:
        values["land_elevation_mean_m"] = 0.0
        values["land_elevation_p95_m"] = 0.0
        values["high_land_fraction_gt2500m"] = 0.0

    terrain_objects = geomorphology or {}
    if terrain_objects:
        group_counts = terrain_objects.get("group_world_feature_counts", {})
        values["geomorphology_world_group_count"] = float(sum(group_counts.values()))
    detail = np.asarray(
        world.get_field("terrain.continental_detail", np.zeros(grid.n)),
        dtype=int,
    )
    if land.any():
        values["orogen_or_plateau_fraction_of_land"] = float(
            np.mean(np.isin(detail[land], [5, 6]))
        )
    else:
        values["orogen_or_plateau_fraction_of_land"] = 0.0

    metrics: dict[str, dict[str, Any]] = {}
    out_of_envelope: list[str] = []
    for name, spec in REFERENCE_ENVELOPES.items():
        lo, hi = spec["range"]
        value = float(values.get(name, 0.0))
        in_range = bool(lo <= value <= hi)
        if not in_range:
            out_of_envelope.append(name)
        metrics[name] = {
            "value": value,
            "reference_min": float(lo),
            "reference_max": float(hi),
            "in_envelope": in_range,
            "basis": spec["basis"],
        }

    acceptance = {
        "screening_only": True,
        "source_is_direct_etopo_raster": False,
        "has_land_ocean_hypsometry": all(
            metrics[key]["in_envelope"]
            for key in (
                "land_fraction",
                "land_elevation_p95_m",
                "shelf_to_abyss_depth_delta_m",
            )
        ),
        "has_ocean_province_distribution": all(
            metrics[key]["in_envelope"]
            for key in (
                "shelf_fraction_of_ocean",
                "slope_rise_fraction_of_ocean",
                "abyss_fraction_of_ocean",
            )
        ),
        "needs_parameter_calibration": bool(out_of_envelope),
    }
    return {
        "schema": "aevum.earth_reference_distribution.v1",
        "status": "needs_calibration" if out_of_envelope else "within_initial_envelope",
        "reference_sources": REFERENCE_SOURCES,
        "metrics": metrics,
        "out_of_envelope": out_of_envelope,
        "acceptance": acceptance,
        "notes": [
            "Initial screening envelopes are deliberately broad.",
            "Direct ETOPO raster sampling should replace these ranges before release-quality claims.",
        ],
    }
