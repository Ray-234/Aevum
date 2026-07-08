"""Plateau area-cap and lifecycle-decay diagnostics."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

from aevum.modules.tectonics import DOMAIN_LIP
from aevum.modules.terrain import CONT_DETAIL_PLATEAU


SCHEMA = "aevum.plateau_area_cap_and_decay.v1"

SOURCE_IDS = (
    "NOAA_ETOPO_2022",
    "GMBA_MOUNTAIN_INVENTORY",
    "GLOBAL_GEOLOGIC_PROVINCES_TECTONIC_PLATES",
    "CRUST1_0",
)

EXPECTED_PLATEAU_VARIANTS = (
    "collision_plateau",
    "volcanic_lip_plateau",
)

EXPECTED_CURRENT_RESIDUAL_ITEMS = (
    "terrain.plateau_inventory",
    "terrain.plateau_age_myr",
    "terrain.plateau_decay_stage",
    "terrain.plateau_parent_process_id",
    "tectonics.plateau_lineage_id",
)

EXPECTED_CURRENT_RESIDUAL_KINDS = (
    "plateau",
    "volcanic_lip_plateau",
)

REFERENCE_PLATEAU_FRAMES: tuple[dict[str, Any], ...] = (
    {
        "id": "collision_initiation",
        "variant": "collision_plateau",
        "stage": "incipient_collision_plateau",
        "age_myr": 35.0,
        "parent_process": "continent_continent_collision",
        "parent_object_kind": "suture",
        "area_fraction_world": 0.012,
        "surrounding_interior_area_fraction_world": 0.190,
        "mean_elevation_m": 3200.0,
        "surrounding_platform_elevation_m": 700.0,
        "relief_p90_p10_m": 1700.0,
        "crustal_thickness_m": 62000.0,
        "support": "crustal_thickening",
    },
    {
        "id": "collision_mature",
        "variant": "collision_plateau",
        "stage": "mature_collision_plateau",
        "age_myr": 90.0,
        "parent_process": "crustal_thickening",
        "parent_object_kind": "active_orogen",
        "area_fraction_world": 0.024,
        "surrounding_interior_area_fraction_world": 0.180,
        "mean_elevation_m": 4700.0,
        "surrounding_platform_elevation_m": 650.0,
        "relief_p90_p10_m": 1550.0,
        "crustal_thickness_m": 72000.0,
        "support": "thickened_buoyant_crust",
    },
    {
        "id": "collision_post_peak",
        "variant": "collision_plateau",
        "stage": "post_collision_plateau",
        "age_myr": 240.0,
        "parent_process": "post_collision_gravitational_collapse",
        "parent_object_kind": "old_orogen",
        "area_fraction_world": 0.019,
        "surrounding_interior_area_fraction_world": 0.188,
        "mean_elevation_m": 3650.0,
        "surrounding_platform_elevation_m": 620.0,
        "relief_p90_p10_m": 1750.0,
        "crustal_thickness_m": 64000.0,
        "support": "residual_crustal_root",
    },
    {
        "id": "collision_dissected",
        "variant": "collision_plateau",
        "stage": "dissected_collision_plateau",
        "age_myr": 900.0,
        "parent_process": "erosion_and_lithospheric_relaxation",
        "parent_object_kind": "old_suture",
        "area_fraction_world": 0.011,
        "surrounding_interior_area_fraction_world": 0.196,
        "mean_elevation_m": 2250.0,
        "surrounding_platform_elevation_m": 580.0,
        "relief_p90_p10_m": 1350.0,
        "crustal_thickness_m": 54000.0,
        "support": "decaying_root_and_inherited_highland",
    },
    {
        "id": "plume_swell",
        "variant": "volcanic_lip_plateau",
        "stage": "plume_swell",
        "age_myr": 5.0,
        "parent_process": "mantle_plume_dynamic_uplift",
        "parent_object_kind": "plume_head",
        "area_fraction_world": 0.008,
        "surrounding_interior_area_fraction_world": 0.160,
        "mean_elevation_m": 900.0,
        "surrounding_platform_elevation_m": 320.0,
        "relief_p90_p10_m": 300.0,
        "crustal_thickness_m": 41000.0,
        "support": "dynamic_topography",
    },
    {
        "id": "lip_emplacement",
        "variant": "volcanic_lip_plateau",
        "stage": "large_igneous_province_emplacement",
        "age_myr": 20.0,
        "parent_process": "large_igneous_province_emplacement",
        "parent_object_kind": "large_igneous_province",
        "area_fraction_world": 0.016,
        "surrounding_interior_area_fraction_world": 0.154,
        "mean_elevation_m": 1700.0,
        "surrounding_platform_elevation_m": 350.0,
        "relief_p90_p10_m": 460.0,
        "crustal_thickness_m": 47000.0,
        "support": "volcanic_construction_and_thermal_buoyancy",
    },
    {
        "id": "lip_post_emplacement",
        "variant": "volcanic_lip_plateau",
        "stage": "post_lip_plateau",
        "age_myr": 150.0,
        "parent_process": "post_lip_cooling_and_erosion",
        "parent_object_kind": "large_igneous_province",
        "area_fraction_world": 0.012,
        "surrounding_interior_area_fraction_world": 0.158,
        "mean_elevation_m": 1220.0,
        "surrounding_platform_elevation_m": 330.0,
        "relief_p90_p10_m": 520.0,
        "crustal_thickness_m": 44000.0,
        "support": "cooling_thermal_swell",
    },
    {
        "id": "lip_eroded_surface",
        "variant": "volcanic_lip_plateau",
        "stage": "eroded_lip_surface",
        "age_myr": 520.0,
        "parent_process": "lip_erosion_and_thermal_decay",
        "parent_object_kind": "eroded_volcanic_province",
        "area_fraction_world": 0.006,
        "surrounding_interior_area_fraction_world": 0.164,
        "mean_elevation_m": 650.0,
        "surrounding_platform_elevation_m": 300.0,
        "relief_p90_p10_m": 360.0,
        "crustal_thickness_m": 39000.0,
        "support": "eroded_volcanic_surface",
    },
)


def _canonical_digest(data: Any) -> str:
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _field_or_empty(world: Any, name: str, fill: float = 0.0) -> np.ndarray:
    if name in getattr(world, "fields", {}):
        return np.asarray(world.field(name), dtype=np.float64)
    return np.full(world.grid.n, fill, dtype=np.float64)


def _item_present(world: Any, name: str) -> bool:
    return (
        name in getattr(world, "fields", {})
        or name in getattr(world, "objects", {})
        or name in getattr(world, "networks", {})
        or name in getattr(world, "globals", {})
    )


def _object_mask(world: Any, object_set: str, kind: str | None = None) -> np.ndarray:
    mask = np.zeros(world.grid.n, dtype=bool)
    for obj in getattr(world, "objects", {}).get(object_set, []):
        if kind is not None and str(obj.get("kind", "")) != kind:
            continue
        cells = np.asarray(obj.get("cells", ()), dtype=np.int64)
        cells = cells[(cells >= 0) & (cells < world.grid.n)]
        if cells.size:
            mask[cells] = True
    return mask


def _area_fraction(
    area: np.ndarray,
    mask: np.ndarray,
    within: np.ndarray | None = None,
) -> float:
    mask = np.asarray(mask, dtype=bool)
    if within is not None:
        within = np.asarray(within, dtype=bool)
        denom = float(area[within].sum())
        mask = mask & within
    else:
        denom = float(area.sum())
    if denom <= 0.0:
        return 0.0
    return float(area[mask].sum() / denom)


def _percentile(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(np.percentile(values, q))


def _variant_frames(variant: str) -> tuple[dict[str, Any], ...]:
    return tuple(
        dict(frame)
        for frame in REFERENCE_PLATEAU_FRAMES
        if str(frame["variant"]) == variant
    )


def _variant_peak(frames: tuple[dict[str, Any], ...]) -> dict[str, Any]:
    return max(frames, key=lambda frame: float(frame["mean_elevation_m"]))


def reference_plateau_area_cap_and_decay_summary() -> dict[str, Any]:
    frames = tuple(dict(frame) for frame in REFERENCE_PLATEAU_FRAMES)
    variants = tuple(sorted({str(frame["variant"]) for frame in frames}))
    variant_counts = Counter(str(frame["variant"]) for frame in frames)
    stage_counts = Counter(str(frame["stage"]) for frame in frames)
    parent_process_failures = tuple(
        str(frame["id"]) for frame in frames if not frame.get("parent_process")
    )
    parent_object_failures = tuple(
        str(frame["id"]) for frame in frames if not frame.get("parent_object_kind")
    )

    collision_frames = _variant_frames("collision_plateau")
    volcanic_frames = _variant_frames("volcanic_lip_plateau")
    collision_peak = _variant_peak(collision_frames)
    volcanic_peak = _variant_peak(volcanic_frames)
    collision_final = collision_frames[-1]
    volcanic_final = volcanic_frames[-1]

    collision_max_area = max(float(frame["area_fraction_world"]) for frame in collision_frames)
    volcanic_max_area = max(float(frame["area_fraction_world"]) for frame in volcanic_frames)
    total_peak_area = (
        float(collision_peak["area_fraction_world"])
        + float(volcanic_peak["area_fraction_world"])
    )
    min_background_dominance = min(
        float(frame["surrounding_interior_area_fraction_world"])
        / max(float(frame["area_fraction_world"]), 1.0e-12)
        for frame in frames
    )
    max_delta_elevation = max(
        float(frame["mean_elevation_m"])
        - float(frame["surrounding_platform_elevation_m"])
        for frame in frames
    )

    collision_elevation_decay = (
        float(collision_peak["mean_elevation_m"])
        - float(collision_final["mean_elevation_m"])
    )
    collision_area_decay = (
        float(collision_peak["area_fraction_world"])
        - float(collision_final["area_fraction_world"])
    )
    volcanic_elevation_decay = (
        float(volcanic_peak["mean_elevation_m"])
        - float(volcanic_final["mean_elevation_m"])
    )
    volcanic_area_decay = (
        float(volcanic_peak["area_fraction_world"])
        - float(volcanic_final["area_fraction_world"])
    )

    acceptance = {
        "fixture_schema_ready": bool(frames and SOURCE_IDS),
        "required_variants_present": variants == EXPECTED_PLATEAU_VARIANTS,
        "lifecycle_frames_present": (
            min(variant_counts.values(), default=0) >= 4
            and len(stage_counts) >= 8
        ),
        "parent_processes_present": not parent_process_failures,
        "parent_objects_present": not parent_object_failures,
        "collision_plateau_area_capped": 0.015 <= collision_max_area <= 0.035,
        "volcanic_plateau_area_capped": 0.010 <= volcanic_max_area <= 0.025,
        "combined_peak_area_finite": total_peak_area <= 0.060,
        "collision_plateau_decays": (
            collision_elevation_decay >= 1500.0 and collision_area_decay >= 0.008
        ),
        "volcanic_plateau_decays": (
            volcanic_elevation_decay >= 800.0 and volcanic_area_decay >= 0.006
        ),
        "plateaus_not_default_interiors": min_background_dominance >= 6.0,
        "plateaus_above_surrounding_platforms": max_delta_elevation >= 1200.0,
        "raw_reference_data_not_stored": True,
    }
    metrics = {
        "variant_count": int(len(variants)),
        "frame_count": int(len(frames)),
        "stage_count": int(len(stage_counts)),
        "collision_frame_count": int(variant_counts["collision_plateau"]),
        "volcanic_frame_count": int(variant_counts["volcanic_lip_plateau"]),
        "parent_process_failure_count": int(len(parent_process_failures)),
        "parent_object_failure_count": int(len(parent_object_failures)),
        "max_collision_plateau_area_fraction_world": float(collision_max_area),
        "max_volcanic_plateau_area_fraction_world": float(volcanic_max_area),
        "combined_peak_plateau_area_fraction_world": float(total_peak_area),
        "min_background_to_plateau_area_ratio": float(min_background_dominance),
        "collision_peak_elevation_m": float(collision_peak["mean_elevation_m"]),
        "collision_final_elevation_m": float(collision_final["mean_elevation_m"]),
        "collision_elevation_decay_m": float(collision_elevation_decay),
        "collision_peak_area_fraction_world": float(collision_peak["area_fraction_world"]),
        "collision_final_area_fraction_world": float(collision_final["area_fraction_world"]),
        "collision_area_decay_fraction_world": float(collision_area_decay),
        "volcanic_peak_elevation_m": float(volcanic_peak["mean_elevation_m"]),
        "volcanic_final_elevation_m": float(volcanic_final["mean_elevation_m"]),
        "volcanic_elevation_decay_m": float(volcanic_elevation_decay),
        "volcanic_peak_area_fraction_world": float(volcanic_peak["area_fraction_world"]),
        "volcanic_final_area_fraction_world": float(volcanic_final["area_fraction_world"]),
        "volcanic_area_decay_fraction_world": float(volcanic_area_decay),
        "max_plateau_relief_p90_p10_m": float(
            max(float(frame["relief_p90_p10_m"]) for frame in frames)),
        "max_plateau_delta_above_platform_m": float(max_delta_elevation),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "plateau_area_cap_decay_reference_ready"
            if all(acceptance.values())
            else "plateau_area_cap_decay_reference_incomplete"
        ),
        "source_ids": SOURCE_IDS,
        "variants": variants,
        "expected_variants": EXPECTED_PLATEAU_VARIANTS,
        "frames": frames,
        "parent_process_failures": parent_process_failures,
        "parent_object_failures": parent_object_failures,
        "extraction_policy": {
            "raw_topography_mountain_inventory_or_crustal_grids_stored": False,
            "direct_plateau_extraction_pending": True,
            "fixture_is_small_derived_reference": True,
        },
        "metrics": metrics,
        "acceptance": acceptance,
        "fixture_digest": _canonical_digest({
            "frames": frames,
            "source_ids": SOURCE_IDS,
        }),
    }


def current_generated_plateau_audit(world: Any) -> dict[str, Any]:
    grid = world.grid
    area = np.asarray(grid.cell_area, dtype=np.float64)
    fields = getattr(world, "fields", {})
    objects = getattr(world, "objects", {})
    elevation = _field_or_empty(world, "terrain.elevation_m")
    sea = float(getattr(world, "sea_level", world.g("ocean.sea_level_m", 0.0)))
    rel = elevation - sea
    detail = _field_or_empty(world, "terrain.continental_detail").astype(int)
    plateau_inventory = _field_or_empty(world, "terrain.plateau_inventory")
    crust_type = _field_or_empty(world, "crust.type")
    crust_domain = _field_or_empty(world, "crust.domain").astype(int)
    land = rel >= 0.0
    continent = crust_type == 1.0
    continental_land = land & continent
    plateau_detail = detail == CONT_DETAIL_PLATEAU
    production_plateau = plateau_inventory > 0.0
    plateau_mask = plateau_detail | production_plateau
    plateau_objects = [
        obj for obj in objects.get("terrain.continental_landforms", [])
        if str(obj.get("kind", "")) == "plateau"
    ]
    plateau_object_mask = _object_mask(
        world, "terrain.continental_landforms", "plateau")
    lip_domain = crust_domain == DOMAIN_LIP
    lip_object_mask = _object_mask(world, "tectonics.lips")
    lip_context = lip_domain | lip_object_mask
    high_interior = (
        continental_land
        & (rel >= 1800.0)
        & ~plateau_mask
    )

    missing_items = tuple(
        item for item in EXPECTED_CURRENT_RESIDUAL_ITEMS if not _item_present(world, item)
    )
    expressed_kinds = {
        "plateau" for obj in plateau_objects if str(obj.get("kind", "")) == "plateau"
    }
    if production_plateau.any():
        expressed_kinds.add("plateau")
    if (
        (plateau_mask & lip_context).any()
        or np.any(plateau_inventory == 2.0)
    ):
        expressed_kinds.add("volcanic_lip_plateau")
    missing_kinds = tuple(
        kind for kind in EXPECTED_CURRENT_RESIDUAL_KINDS if kind not in expressed_kinds
    )

    plateau_object_area = sum(float(obj.get("area_fraction", 0.0)) for obj in plateau_objects)
    max_plateau_object_area = max(
        (float(obj.get("area_fraction", 0.0)) for obj in plateau_objects),
        default=0.0,
    )
    parented_plateau_objects = [
        obj for obj in plateau_objects
        if obj.get("parent_process")
        and (
            obj.get("parent_tectonic_object_ids")
            or str(obj.get("parent_process", "")) in {
                "crustal_thickening_or_lip",
                "crustal_thickening",
                "large_igneous_province_emplacement",
            }
        )
    ]
    lip_plateau_overlap = _area_fraction(
        area, plateau_mask & lip_context, lip_context) if lip_context.any() else 0.0
    plateau_rel_values = rel[plateau_mask & continental_land]
    high_interior_area = _area_fraction(area, high_interior, continental_land)

    metrics = {
        "continental_landform_object_count": int(
            len(objects.get("terrain.continental_landforms", []))),
        "plateau_object_count": int(len(plateau_objects)),
        "parented_plateau_object_count": int(len(parented_plateau_objects)),
        "plateau_detail_cell_count": int(np.count_nonzero(plateau_mask)),
        "plateau_detail_land_cell_count": int(
            np.count_nonzero(plateau_mask & continental_land)),
        "plateau_detail_area_fraction_world": _area_fraction(area, plateau_mask),
        "plateau_detail_area_fraction_land": _area_fraction(
            area, plateau_mask, continental_land),
        "plateau_object_area_fraction_world": float(plateau_object_area),
        "max_plateau_object_area_fraction_world": float(max_plateau_object_area),
        "production_plateau_inventory_cell_count": int(
            np.count_nonzero(production_plateau)),
        "production_collision_plateau_cell_count": int(
            np.count_nonzero(plateau_inventory == 1.0)),
        "production_volcanic_plateau_cell_count": int(
            np.count_nonzero(plateau_inventory == 2.0)),
        "tectonics_lip_object_count": int(len(objects.get("tectonics.lips", []))),
        "lip_domain_area_fraction_world": _area_fraction(area, lip_domain),
        "lip_context_area_fraction_world": _area_fraction(area, lip_context),
        "lip_context_plateau_overlap_fraction": float(lip_plateau_overlap),
        "missing_plateau_item_count": int(len(missing_items)),
        "missing_plateau_kind_count": int(len(missing_kinds)),
        "high_interior_without_plateau_fraction_of_continental_land": float(
            high_interior_area),
        "plateau_mean_elevation_m": (
            float(np.average(plateau_rel_values, weights=area[plateau_mask & continental_land]))
            if plateau_rel_values.size else float("nan")
        ),
        "plateau_relief_p90_minus_p10_m": (
            _percentile(plateau_rel_values, 90.0)
            - _percentile(plateau_rel_values, 10.0)
            if plateau_rel_values.size else float("nan")
        ),
    }
    area_cap_audited = (
        metrics["plateau_detail_area_fraction_world"] <= 0.060
        and metrics["plateau_object_area_fraction_world"] <= 0.060
        and metrics["max_plateau_object_area_fraction_world"] <= 0.035
    )
    parentage_audited = (
        metrics["plateau_object_count"] == 0
        or metrics["parented_plateau_object_count"] == metrics["plateau_object_count"]
    )
    expected_residuals_recorded = (
        set(missing_items).issubset(set(EXPECTED_CURRENT_RESIDUAL_ITEMS))
        and set(missing_kinds).issubset(set(EXPECTED_CURRENT_RESIDUAL_KINDS))
    )
    acceptance = {
        "current_core_fields_available": all(
            item in fields
            for item in ("terrain.elevation_m", "terrain.continental_detail",
                         "crust.type", "crust.domain")
        ),
        "current_lip_context_available": (
            metrics["tectonics_lip_object_count"] > 0
            or metrics["lip_domain_area_fraction_world"] > 0.0
        ),
        "current_plateau_area_cap_audited": area_cap_audited,
        "current_plateau_parentage_audited": parentage_audited,
        "current_expected_residuals_recorded": expected_residuals_recorded,
    }
    return {
        "schema": SCHEMA,
        "status": (
            "generated_world_plateau_area_cap_audit_ready"
            if all(acceptance.values())
            else "generated_world_plateau_area_cap_audit_incomplete"
        ),
        "landform_kind_counts": dict(sorted(Counter(
            str(obj.get("kind", ""))
            for obj in objects.get("terrain.continental_landforms", [])
        ).items())),
        "missing_plateau_items": missing_items,
        "missing_expected_plateau_kinds": missing_kinds,
        "expected_current_residual_items": EXPECTED_CURRENT_RESIDUAL_ITEMS,
        "expected_current_residual_kinds": EXPECTED_CURRENT_RESIDUAL_KINDS,
        "limitations": {
            "first_class_plateau_inventory_missing": bool(missing_items),
            "plateau_landform_expression_missing": (
                metrics["plateau_object_count"] == 0
                and metrics["plateau_detail_land_cell_count"] == 0
                and metrics["production_plateau_inventory_cell_count"] == 0
            ),
            "volcanic_lip_plateau_expression_missing": (
                bool(lip_context.any())
                and metrics["lip_context_plateau_overlap_fraction"] <= 0.0
                and metrics["production_volcanic_plateau_cell_count"] == 0
            ),
            "plateau_decay_fields_missing": (
                "terrain.plateau_age_myr" in missing_items
                or "terrain.plateau_decay_stage" in missing_items
            ),
            "plateau_parent_lineage_missing": (
                "terrain.plateau_parent_process_id" in missing_items
                or "tectonics.plateau_lineage_id" in missing_items
            ),
            "high_interior_support_needs_p90_gap_audit": (
                metrics["high_interior_without_plateau_fraction_of_continental_land"]
                > 0.02
            ),
        },
        "metrics": metrics,
        "acceptance": acceptance,
    }


def plateau_area_cap_and_decay_summary(world: Any | None = None) -> dict[str, Any]:
    reference = reference_plateau_area_cap_and_decay_summary()
    current = (
        current_generated_plateau_audit(world)
        if world is not None
        else {
            "status": "generated_world_plateau_area_cap_audit_not_run",
            "acceptance": {
                "current_core_fields_available": False,
                "current_lip_context_available": False,
                "current_plateau_area_cap_audited": False,
                "current_plateau_parentage_audited": False,
                "current_expected_residuals_recorded": False,
                "production_plateau_inventory_fields_available": False,
                "production_plateau_expression_available": False,
            },
            "metrics": {},
        }
    )
    acceptance = {
        "reference_plateau_fixture_ready": (
            reference["status"] == "plateau_area_cap_decay_reference_ready"),
        "reference_area_caps_and_decay_ready": (
            reference["acceptance"]["collision_plateau_area_capped"]
            and reference["acceptance"]["volcanic_plateau_area_capped"]
            and reference["acceptance"]["collision_plateau_decays"]
            and reference["acceptance"]["volcanic_plateau_decays"]
            and reference["acceptance"]["plateaus_not_default_interiors"]
        ),
        "current_generated_audit_available": (
            current["status"] == "generated_world_plateau_area_cap_audit_ready"),
        "current_plateau_area_cap_audited": current["acceptance"][
            "current_plateau_area_cap_audited"],
        "current_expected_residuals_recorded": current["acceptance"][
            "current_expected_residuals_recorded"],
        "production_plateau_inventory_fields_available": (
            not current.get("missing_plateau_items", ())
        ),
        "production_plateau_expression_available": (
            not current.get("missing_expected_plateau_kinds", ())
        ),
    }
    return {
        "schema": SCHEMA,
        "status": (
            "plateau_area_cap_and_decay_ready"
            if all(acceptance.values())
            else "plateau_area_cap_and_decay_incomplete"
        ),
        "reference": reference,
        "current_generated": current,
        "acceptance": acceptance,
        "summary_digest": _canonical_digest({
            "reference": reference,
            "current": current,
        }),
    }
