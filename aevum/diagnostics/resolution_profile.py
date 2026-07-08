"""Resolution-ladder profiling for high-resolution geomorphology reviews.

This diagnostic intentionally stays lighter than the P12 release gate.  It
records build/run/compile/diagnostic costs and key morphology metrics at one or
more cell counts so 24000/72000-cell audits can be planned from evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import importlib.metadata
import importlib.util
import json
import platform
import time
from pathlib import Path
from typing import Any, Callable, TypeVar

import numpy as np

from aevum.compiler.map_compiler import MapCompiler
from aevum.engine import Engine
from aevum.spec.presets import get_preset


T = TypeVar("T")


@dataclass(frozen=True)
class ResolutionProfileConfig:
    preset: str = "earthlike"
    cells: tuple[int, ...] = (900, 2500)
    t_end_myr: float | None = None
    frames: int = 4
    hex_width: int = 64
    hex_height: int = 32
    starts: int = 4
    compile_map: bool = True
    tectonic_diagnostics: bool = True
    geomorphology_coverage: bool = False
    render_assets: bool = False
    progress: bool = False
    projection_cells: tuple[int, ...] = (8000, 24000, 72000)


def run_resolution_profile(
    config: ResolutionProfileConfig,
    outdir: Path | None = None,
) -> dict[str, Any]:
    """Run the same profile workload over a cell-count ladder."""
    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)

    entries = [_run_one(cells, config, outdir) for cells in config.cells]
    summary = {
        "schema": "aevum.resolution_profile.v1",
        "config": {
            "preset": config.preset,
            "cells": [int(c) for c in config.cells],
            "t_end_myr": (
                None if config.t_end_myr is None else float(config.t_end_myr)
            ),
            "frames": int(config.frames),
            "hex_width": int(config.hex_width),
            "hex_height": int(config.hex_height),
            "starts": int(config.starts),
            "compile_map": bool(config.compile_map),
            "tectonic_diagnostics": bool(config.tectonic_diagnostics),
            "geomorphology_coverage": bool(config.geomorphology_coverage),
            "render_assets": bool(config.render_assets),
            "projection_cells": [int(c) for c in config.projection_cells],
        },
        "environment": _environment_summary(),
        "entries": entries,
    }
    summary["scaling"] = _scaling_summary(entries)
    summary["high_resolution"] = _high_resolution_summary(
        entries,
        summary["scaling"],
        config.projection_cells,
        summary["environment"],
    )

    if outdir is not None:
        (outdir / "resolution_profile_summary.json").write_text(
            json.dumps(summary, indent=2, default=_json_default) + "\n"
        )
    return summary


def _run_one(
    cells: int,
    config: ResolutionProfileConfig,
    outdir: Path | None,
) -> dict[str, Any]:
    from aevum import validation

    stage_seconds: dict[str, float] = {}

    spec = get_preset(config.preset)
    spec.grid_cells = int(cells)
    if config.t_end_myr is not None:
        spec.t_end_myr = float(config.t_end_myr)

    eng = _time_stage(stage_seconds, "build", lambda: Engine.build(spec))
    _time_stage(
        stage_seconds,
        "run",
        lambda: eng.run(n_frames=config.frames, progress=config.progress),
    )

    cm = None
    compiler_metrics: dict[str, Any] = {}
    if config.compile_map:
        cm = _time_stage(
            stage_seconds,
            "compile",
            lambda: MapCompiler(eng.world, eng.archive).compile(
                width=config.hex_width,
                height=config.hex_height,
                n_starts=config.starts,
            ),
        )
        compiler_metrics = _time_stage(
            stage_seconds,
            "compiler_diagnostics",
            lambda: validation.compiler_consistency_metrics(eng, cm),
        )

    tectonics: dict[str, Any] | None = None
    if config.tectonic_diagnostics:
        tectonics = _time_stage(
            stage_seconds,
            "tectonic_diagnostics",
            lambda: validation.tectonic_diagnostics(eng),
        )

    coverage: dict[str, Any] | None = None
    if config.geomorphology_coverage:
        coverage = _time_stage(
            stage_seconds,
            "geomorphology_coverage",
            lambda: validation.earth_geomorphology_coverage_metrics(eng, cm),
        )

    rendered: list[str] = []
    if config.render_assets and outdir is not None:
        rendered = _time_stage(
            stage_seconds,
            "render",
            lambda: _render_assets(eng, cm, outdir / f"{config.preset}_{cells}cells"),
        )

    entry: dict[str, Any] = {
        "preset": config.preset,
        "seed": int(spec.seed),
        "cells": int(cells),
        "t_end_myr": float(spec.t_end_myr),
        "time_myr": float(eng.world.time_myr),
        "stage_seconds": stage_seconds,
        "total_seconds": float(sum(stage_seconds.values())),
        "scheduler": _scheduler_summary(eng),
        "world": _world_summary(eng),
        "morphology": _morphology_summary(tectonics),
        "crust": _crust_summary(tectonics),
        "ocean_geography": _ocean_summary(tectonics),
        "terrain": _terrain_summary(eng.world),
        "compiler": _compiler_summary(cm, compiler_metrics),
        "geomorphology_coverage": _coverage_summary(coverage),
        "warnings": list(tectonics.get("warnings", []) if tectonics else []),
        "rendered_assets": rendered,
    }
    return entry


def _time_stage(
    stage_seconds: dict[str, float],
    name: str,
    fn: Callable[[], T],
) -> T:
    t0 = time.perf_counter()
    result = fn()
    stage_seconds[name] = float(time.perf_counter() - t0)
    return result


def _environment_summary() -> dict[str, Any]:
    try:
        import scipy
        scipy_version = scipy.__version__
    except Exception:  # pragma: no cover - diagnostic fallback only
        scipy_version = None
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "numpy": np.__version__,
        "scipy": scipy_version,
        "active_array_backend": "numpy",
        "optional_acceleration": {
            "numba": _optional_package("numba", "numba", ("numba",)),
            "cupy": _optional_package(
                "cupy", "cupy", ("cupy", "cupy-cuda11x", "cupy-cuda12x")
            ),
            "jax": _optional_package("jax", "jax", ("jax",)),
            "torch": _optional_package("torch", "torch", ("torch",)),
        },
    }


def _optional_package(
    name: str,
    module_name: str,
    distribution_names: tuple[str, ...],
) -> dict[str, Any]:
    module_available = importlib.util.find_spec(module_name) is not None
    version = None
    distribution = None
    for dist in distribution_names:
        try:
            version = importlib.metadata.version(dist)
            distribution = dist
            break
        except importlib.metadata.PackageNotFoundError:
            continue
    return {
        "available": bool(module_available),
        "distribution": distribution,
        "version": version,
        "role": _acceleration_role(name),
    }


def _acceleration_role(name: str) -> str:
    roles = {
        "numba": (
            "candidate for graph BFS, width/component passes, articulation, "
            "and priority-flood style loops"
        ),
        "cupy": "candidate for dense field algebra and smoothing kernels",
        "jax": "candidate for dense field algebra behind an explicit backend boundary",
        "torch": "diagnostic only unless a tensor backend is introduced explicitly",
    }
    return roles.get(name, "optional acceleration package")


def _scheduler_summary(eng: Engine) -> dict[str, Any]:
    run_counts = {sm.module.name: int(sm.runs) for sm in eng.scheduler.modules}
    history_run_counts: dict[str, int] = {}
    for rec in eng.scheduler.history:
        for name in rec.ran:
            history_run_counts[name] = history_run_counts.get(name, 0) + 1
    return {
        "macro_steps": int(len(eng.scheduler.history)),
        "archive_frames": int(len(eng.archive.frames)),
        "events": int(len(eng.bus.events)),
        "module_run_counts": run_counts,
        "history_module_run_counts": history_run_counts,
    }


def _world_summary(eng: Engine) -> dict[str, Any]:
    w = eng.world
    return {
        "land_fraction": float(w.land_fraction()),
        "sea_level_m": float(w.sea_level),
        "surface_area_m2": float(w.spec.surface_area_m2),
        "regime": w.spec.initial_tectonic_regime.value,
    }


def _morphology_summary(tectonics: dict[str, Any] | None) -> dict[str, Any]:
    if not tectonics:
        return {}
    land = tectonics["morphology"]["exposed_land"]
    cont = tectonics["morphology"]["continental_crust"]
    coupling = tectonics["morphology"]["crust_land_coupling"]
    return {
        "land_component_count": int(land["component_count"]),
        "largest_land_component_fraction": float(
            land["largest_component_area_fraction_of_mask"]
        ),
        "land_width_p50_steps": float(land["width_p50_steps"]),
        "land_width_p90_steps": float(land["width_p90_steps"]),
        "land_ribbon_fraction_gt_0_5": float(land["ribbon_area_fraction_gt_0_5"]),
        "land_narrow_necks_per_1000": float(
            land["narrow_neck_cells_per_1000_mask_cells"]
        ),
        "land_coastline_complexity_largest": float(
            land["coastline_complexity_largest_component"]
        ),
        "continental_component_count": int(cont["component_count"]),
        "continental_width_p50_steps": float(cont["width_p50_steps"]),
        "continental_width_p90_steps": float(cont["width_p90_steps"]),
        "continental_ribbon_fraction_gt_0_5": float(
            cont["ribbon_area_fraction_gt_0_5"]
        ),
        "continental_narrow_fraction_le2_width": float(
            cont["narrow_fraction_le2_width"]
        ),
        "exposed_continental_land_fraction": float(
            coupling["exposed_continental_land_fraction_of_land"]
        ),
        "high_oceanic_land_fraction_gt1500m": float(
            coupling["high_oceanic_land_fraction_gt1500m_of_land"]
        ),
    }


def _crust_summary(tectonics: dict[str, Any] | None) -> dict[str, Any]:
    if not tectonics:
        return {}
    crust = tectonics["crust_distribution"]
    return {
        "continental_area_fraction": float(crust["continental_area_fraction"]),
        "stable_craton_fraction_gt075": float(crust["stable_craton_fraction_gt075"]),
        "ancient_continental_fraction_gt2500_myr": float(
            crust["ancient_continental_fraction_gt2500_myr"]
        ),
        "oceanic_age_p95_myr": float(crust["oceanic_age_p95_myr"]),
        "oceanic_old_fraction_gt300_myr": float(crust["oceanic_old_fraction_gt300_myr"]),
    }


def _ocean_summary(tectonics: dict[str, Any] | None) -> dict[str, Any]:
    if not tectonics:
        return {}
    ocean = tectonics["ocean_geography"]
    return {
        "basin_count": int(ocean["ocean_basin_count"]),
        "shelf_fraction": float(ocean["shelf_fraction_of_ocean"]),
        "slope_rise_fraction": float(ocean["slope_rise_fraction_of_ocean"]),
        "abyss_fraction": float(ocean["abyss_fraction_of_ocean"]),
        "ridge_fraction": float(ocean["ridge_fraction_of_ocean"]),
        "trench_fraction": float(ocean["trench_fraction_of_ocean"]),
        "nearshore_depth_p75_m": float(ocean["nearshore_depth_p75_m"]),
        "shelf_to_abyss_depth_delta_m": float(ocean["shelf_to_abyss_depth_delta_m"]),
        "nearshore_superdeep_fraction_gt2500m": float(
            ocean["nearshore_superdeep_fraction_gt2500m"]
        ),
        "far_ocean_shallow_fraction_lt1500m": float(
            ocean["far_ocean_shallow_fraction_lt1500m"]
        ),
    }


def _terrain_summary(world) -> dict[str, Any]:
    elev = np.asarray(world.get_field("terrain.elevation_m"), dtype=np.float64)
    rel = elev - float(world.sea_level)
    land = rel >= 0.0
    if not land.any():
        return {
            "land_elevation_mean_m": 0.0,
            "land_elevation_p75_m": 0.0,
            "land_elevation_p95_m": 0.0,
            "land_elevation_max_m": 0.0,
        }
    land_rel = rel[land]
    return {
        "land_elevation_mean_m": float(np.mean(land_rel)),
        "land_elevation_p75_m": float(np.percentile(land_rel, 75)),
        "land_elevation_p95_m": float(np.percentile(land_rel, 95)),
        "land_elevation_max_m": float(np.max(land_rel)),
    }


def _compiler_summary(cm, compiler_metrics: dict[str, Any]) -> dict[str, Any]:
    if cm is None:
        return {"compiled": False}
    return {
        "compiled": True,
        "width": int(cm.meta["width"]),
        "height": int(cm.meta["height"]),
        "starts": int(len(cm.starts)),
        "compiled_land_fraction": float(
            compiler_metrics.get("compiled_land_fraction", cm.meta["land_fraction"])
        ),
        "source_land_fraction_mean": float(
            compiler_metrics.get("source_land_fraction_mean", 0.0)
        ),
        "land_fraction_abs_delta_from_source": float(
            compiler_metrics.get("land_fraction_abs_delta_from_source", 0.0)
        ),
        "terrain_elevation_sign_mismatch_fraction": float(
            compiler_metrics.get("terrain_elevation_sign_mismatch_fraction", 0.0)
        ),
        "passed_envelope": bool(compiler_metrics.get("passed_envelope", False)),
        "local_yield_cv": float(cm.fairness.get("local_yield_cv", 0.0)),
    }


def _coverage_summary(coverage: dict[str, Any] | None) -> dict[str, Any]:
    if not coverage:
        return {"computed": False}
    return {
        "computed": True,
        "covered_feature_count": int(coverage["covered_feature_count"]),
        "partial_feature_count": int(coverage["partial_feature_count"]),
        "weak_feature_count": int(coverage["weak_feature_count"]),
        "missing_feature_count": int(coverage["missing_feature_count"]),
        "parentless_major_landform_fraction": float(
            coverage["parentless_major_landform_fraction"]
        ),
        "group_world_feature_counts": coverage["group_world_feature_counts"],
        "hard_failures": list(coverage.get("hard_failures", [])),
        "warnings": list(coverage.get("warnings", [])),
    }


def _scaling_summary(entries: list[dict[str, Any]]) -> dict[str, Any]:
    if len(entries) < 2:
        return {"computed": False}
    by_stage: dict[str, list[dict[str, float]]] = {}
    prev = entries[0]
    for cur in entries[1:]:
        cell_ratio = float(cur["cells"]) / max(float(prev["cells"]), 1.0)
        stages = set(prev["stage_seconds"]) | set(cur["stage_seconds"])
        for stage in sorted(stages):
            a = float(prev["stage_seconds"].get(stage, 0.0))
            b = float(cur["stage_seconds"].get(stage, 0.0))
            if a <= 0.0 or b <= 0.0 or cell_ratio <= 1.0:
                exponent = None
            else:
                exponent = float(np.log(b / a) / np.log(cell_ratio))
            by_stage.setdefault(stage, []).append({
                "from_cells": int(prev["cells"]),
                "to_cells": int(cur["cells"]),
                "seconds_ratio": float(b / a) if a > 0.0 else 0.0,
                "cell_ratio": cell_ratio,
                "scaling_exponent": exponent,
            })
        prev = cur
    return {"computed": True, "stage_scaling": by_stage}


def _high_resolution_summary(
    entries: list[dict[str, Any]],
    scaling: dict[str, Any],
    projection_cells: tuple[int, ...],
    environment: dict[str, Any],
) -> dict[str, Any]:
    if not entries:
        return {"computed": False}

    largest = max(entries, key=lambda e: int(e["cells"]))
    largest_cells = int(largest["cells"])
    stage_seconds = {
        str(stage): float(seconds)
        for stage, seconds in largest.get("stage_seconds", {}).items()
    }
    bottlenecks = [
        {"stage": stage, "seconds": seconds}
        for stage, seconds in sorted(
            stage_seconds.items(), key=lambda item: item[1], reverse=True
        )
    ]

    projections = []
    for target in sorted({int(c) for c in projection_cells if int(c) > largest_cells}):
        stage_estimates: dict[str, float] = {}
        ratio = float(target) / max(float(largest_cells), 1.0)
        for stage, seconds in stage_seconds.items():
            exponent = _latest_stage_exponent(scaling, stage)
            if exponent is None:
                exponent = _default_stage_exponent(stage)
            exponent = float(np.clip(exponent, 0.5, 2.5))
            stage_estimates[stage] = float(seconds * ratio ** exponent)
        projections.append({
            "target_cells": target,
            "from_cells": largest_cells,
            "estimated_total_seconds": float(sum(stage_estimates.values())),
            "stage_estimated_seconds": stage_estimates,
            "estimation_note": (
                "coarse extrapolation from the largest measured profile entry; "
                "validate with the next resolution tier before relying on it"
            ),
        })

    installed = [
        name
        for name, meta in environment.get("optional_acceleration", {}).items()
        if meta.get("available")
    ]
    readiness = "preflight_only"
    if largest_cells >= 24000:
        readiness = "medium_deployment_profiled"
    elif largest_cells >= 8000:
        readiness = "routine_review_profiled"

    return {
        "computed": True,
        "readiness": readiness,
        "measured_max_cells": largest_cells,
        "resolution_tiers": [
            {
                "cells": "900-2500",
                "role": "CI/microbenchmark scale; mechanism tests and fast regressions",
            },
            {
                "cells": "8000",
                "role": "routine Earth-like map review and major topology regression",
            },
            {
                "cells": "24000",
                "role": "medium deployment review for coastlines, arcs, shelves, and deltas",
            },
            {
                "cells": "72000",
                "role": "occasional offline audit for isthmuses, small mountains, straits, and delta/fan geometry",
            },
        ],
        "bottleneck_stages_at_max_resolution": bottlenecks,
        "projection_estimates": projections,
        "active_array_backend": environment.get("active_array_backend", "numpy"),
        "installed_optional_acceleration": installed,
        "acceleration_policy": (
            "Keep NumPy/SciPy as the reference path; optimize measured CPU "
            "hotspots before adding optional Numba/CuPy/JAX backends."
        ),
        "optimization_order": [
            "cache common neighbour degree, width, and connected-component results",
            "replace repeated Python neighbour loops with edge-array reductions where equivalent",
            "profile terrain smoothing and hydrology routing before GPU work",
            "add Numba only for graph-heavy loops with parity tests",
            "add GPU array backends only behind explicit backend boundaries for dense kernels",
        ],
    }


def _latest_stage_exponent(scaling: dict[str, Any], stage: str) -> float | None:
    if not scaling.get("computed"):
        return None
    series = scaling.get("stage_scaling", {}).get(stage, [])
    for item in reversed(series):
        exponent = item.get("scaling_exponent")
        if exponent is not None and np.isfinite(exponent):
            return float(exponent)
    return None


def _default_stage_exponent(stage: str) -> float:
    defaults = {
        "build": 1.15,
        "run": 1.25,
        "compile": 1.05,
        "compiler_diagnostics": 1.05,
        "tectonic_diagnostics": 1.15,
        "geomorphology_coverage": 1.15,
        "render": 1.05,
    }
    return defaults.get(stage, 1.20)


def _render_assets(eng: Engine, cm, outdir: Path) -> list[str]:
    from aevum import render

    outdir.mkdir(parents=True, exist_ok=True)
    paths = list(render.render_world(eng.world, outdir))
    if cm is not None:
        paths.append(render.render_hexmap(cm, outdir))
        paths.append(render.render_compiler_consistency(cm, outdir))
    return [str(path) for path in paths]


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)
