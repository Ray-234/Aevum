"""Integrated real-Earth morphology promotion audit helpers."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

from aevum.diagnostics.current_world_morphology_gap_inventory import (
    REQUIRED_REVIEW_ASSETS,
)


SCHEMA = "aevum.integrated_real_earth_morphology_promotion_audit.v1"

STAGE_SUITES = tuple(f"P{idx}" for idx in range(76, 91))

PROMOTION_HIERARCHY = (
    "process_object_correctness",
    "archive_lifecycle_continuity",
    "province_crust_sediment_drainage_ordering",
    "hypsometry_bathymetry_planform_envelopes",
    "compiler_render_agreement",
    "png_contact_sheet_visual_review",
)


def _resolve_path(root: Path, path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return root / path


def latest_stage_summary_matrix(
    root: Path,
    *,
    current_summaries: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Read the latest archived P76-P90 benchmark summaries.

    The matrix is evidence for the integrated audit; it deliberately does not
    rerun every earlier generated-world benchmark.
    """
    current_summaries = current_summaries or {}
    rows: dict[str, dict[str, Any]] = {}
    for suite in STAGE_SUITES:
        if suite in current_summaries:
            summary = current_summaries[suite]
            benches = list(summary.get("benchmarks", ()))
            rows[suite] = {
                "suite": suite,
                "source": "current",
                "summary_path": "",
                "schema": str(summary.get("schema", "")),
                "status": str(summary.get("status", "")),
                "benchmark_count": int(len(benches)),
                "passed_benchmark_count": int(
                    sum(1 for bench in benches if bench.get("passed"))),
                "all_microbenchmarks_pass": bool(
                    summary.get("acceptance", {}).get(
                        "all_microbenchmarks_pass",
                        summary.get("status") == "pass",
                    )
                ),
            }
            continue

        pattern = f"out_bench_{suite.lower()}_*/tectonics_bench_summary.json"
        candidates = sorted(root.glob(pattern))
        selected_summary: dict[str, Any] | None = None
        selected_path: Path | None = None
        load_error = ""
        for path in reversed(candidates):
            try:
                loaded = json.loads(path.read_text())
            except Exception as exc:
                load_error = str(exc)
                continue
            if str(loaded.get("suite", "")).upper() == suite:
                selected_summary = loaded
                selected_path = path
                break
        if selected_summary is None:
            rows[suite] = {
                "suite": suite,
                "source": "missing",
                "summary_path": "",
                "schema": "",
                "status": "missing",
                "benchmark_count": 0,
                "passed_benchmark_count": 0,
                "all_microbenchmarks_pass": False,
                "load_error": load_error,
            }
            continue
        benches = list(selected_summary.get("benchmarks", ()))
        rows[suite] = {
            "suite": suite,
            "source": "archived",
            "summary_path": str(selected_path),
            "schema": str(selected_summary.get("schema", "")),
            "status": str(selected_summary.get("status", "")),
            "benchmark_count": int(len(benches)),
            "passed_benchmark_count": int(
                sum(1 for bench in benches if bench.get("passed"))),
            "all_microbenchmarks_pass": bool(
                selected_summary.get("acceptance", {}).get(
                    "all_microbenchmarks_pass",
                    selected_summary.get("status") == "pass",
                )
            ),
        }

    missing = tuple(suite for suite, row in rows.items()
                    if row["status"] == "missing")
    failing = tuple(suite for suite, row in rows.items()
                    if row["status"] != "pass"
                    or not row["all_microbenchmarks_pass"])
    return {
        "schema": "aevum.p91_stage_summary_matrix.v1",
        "suite_count": int(len(rows)),
        "pass_count": int(sum(
            row["status"] == "pass" and row["all_microbenchmarks_pass"]
            for row in rows.values()
        )),
        "missing_count": int(len(missing)),
        "failing_count": int(len(failing)),
        "missing_suites": missing,
        "failing_suites": failing,
        "rows": rows,
    }


def highres_review_asset_metrics(
    root: Path,
    p69_summary: dict[str, Any] | None,
) -> dict[str, Any]:
    if not p69_summary:
        return {
            "available": False,
            "p69_status": "",
            "p69_cells": 0,
            "member_count": 0,
            "required_asset_count": int(len(REQUIRED_REVIEW_ASSETS)),
            "expected_required_asset_files": 0,
            "actual_required_asset_files": 0,
            "member_asset_dirs_present": 0,
            "complete_member_count": 0,
            "missing_assets_by_member": {},
            "contact_sheet_paths": (),
            "contact_sheet_present_count": 0,
            "all_required_assets_present": False,
        }

    metrics = p69_summary.get("benchmarks", [{}])[0].get("metrics", {})
    rows = list(metrics.get("per_member_metrics", ()))
    missing_by_member: dict[str, tuple[str, ...]] = {}
    actual = 0
    asset_dirs = 0
    complete = 0
    contact_sheets: set[str] = set()
    for row in rows:
        member = str(row.get("member_name", "unknown_member"))
        asset_dir_text = str(row.get("asset_dir", ""))
        missing: list[str] = []
        if asset_dir_text and _resolve_path(root, asset_dir_text).exists():
            asset_dirs += 1
        for asset in REQUIRED_REVIEW_ASSETS:
            path = _resolve_path(root, str(Path(asset_dir_text) / asset))
            if path.exists():
                actual += 1
            else:
                missing.append(asset)
        if not missing:
            complete += 1
        missing_by_member[member] = tuple(missing)
        contact = str(row.get("contact_sheet", ""))
        if contact:
            contact_sheets.add(contact)

    contact_paths = tuple(sorted(contact_sheets))
    contact_present = int(
        sum(_resolve_path(root, path).exists() for path in contact_paths))
    expected = int(len(rows) * len(REQUIRED_REVIEW_ASSETS))
    return {
        "available": bool(p69_summary),
        "p69_status": str(p69_summary.get("status", "")),
        "p69_cells": int(metrics.get("cells", 0)),
        "p69_highres_ready_all_members": bool(
            p69_summary.get("acceptance", {}).get("highres_ready_all_members", False)),
        "p69_promotion_ready": bool(
            p69_summary.get("acceptance", {}).get("promotion_ready", False)),
        "p69_promotion_blockers": tuple(metrics.get("promotion_blockers", ())),
        "p69_earth_reference_out_of_envelope_max": int(
            metrics.get("earth_reference_out_of_envelope_max", 0)),
        "member_count": int(len(rows)),
        "required_asset_count": int(len(REQUIRED_REVIEW_ASSETS)),
        "expected_required_asset_files": expected,
        "actual_required_asset_files": int(actual),
        "member_asset_dirs_present": int(asset_dirs),
        "complete_member_count": int(complete),
        "missing_assets_by_member": missing_by_member,
        "contact_sheet_paths": contact_paths,
        "contact_sheet_present_count": contact_present,
        "all_required_assets_present": bool(
            expected > 0
            and actual == expected
            and asset_dirs == len(rows)
            and contact_present >= 1
        ),
    }


def ci_world_review_metrics(rows: list[dict[str, Any]],
                            contact_sheet_path: Path | None) -> dict[str, Any]:
    required = len(REQUIRED_REVIEW_ASSETS)
    complete = int(sum(row.get("asset_set_complete", False) for row in rows))
    compiler_passed = int(sum(
        row.get("compiler_passed_envelope", False) for row in rows))
    inventory_ready = int(sum(
        row.get("inventory_status") == "current_world_morphology_gap_inventory_ready"
        for row in rows
    ))
    missing_asset_total = int(sum(
        len(row.get("missing_review_assets", ())) for row in rows))
    contact_present = bool(contact_sheet_path and contact_sheet_path.exists())
    return {
        "schema": "aevum.p91_ci_world_review.v1",
        "world_count": int(len(rows)),
        "cells": tuple(int(row.get("cells", 0)) for row in rows),
        "required_asset_count": int(required),
        "expected_required_asset_files": int(len(rows) * required),
        "actual_required_asset_files": int(sum(
            row.get("required_asset_present_count", 0) for row in rows)),
        "asset_set_complete_count": complete,
        "compiler_passed_count": compiler_passed,
        "inventory_ready_count": inventory_ready,
        "missing_required_asset_count": missing_asset_total,
        "contact_sheet_path": "" if contact_sheet_path is None else str(contact_sheet_path),
        "contact_sheet_present": contact_present,
        "all_required_assets_present": bool(complete == len(rows) and missing_asset_total == 0),
        "all_compilers_passed": bool(compiler_passed == len(rows)),
        "all_gap_inventories_ready": bool(inventory_ready == len(rows)),
    }


def render_ci_world_contact_sheet(rows: list[dict[str, Any]], path: Path) -> Path | None:
    if not rows:
        return None

    from aevum import render

    path.parent.mkdir(parents=True, exist_ok=True)
    panel_specs = (
        "elevation",
        "terrain_provinces",
        "continental_detail",
        "ocean_depth",
        "crust_age",
        "hex_terrain",
    )
    fig, axes = plt.subplots(
        len(rows),
        len(panel_specs),
        figsize=(3.0 * len(panel_specs), 2.35 * len(rows)),
        squeeze=False,
        constrained_layout=True,
    )
    for row_i, row in enumerate(rows):
        world = row["world"]
        compiled = row["compiled"]
        grid = world.grid
        sea = world.sea_level
        fields = {
            "elevation": np.asarray(world.get_field("terrain.elevation_m") - sea),
            "terrain_provinces": np.asarray(world.get_field("terrain.province", 0.0)),
            "continental_detail": np.asarray(
                world.get_field("terrain.continental_detail", 0.0)),
            "ocean_depth": np.asarray(world.get_field("ocean.depth_province", 0.0)),
            "crust_age": np.asarray(world.get_field("crust.age_myr", 0.0)),
        }
        for col_i, name in enumerate(panel_specs):
            ax = axes[row_i, col_i]
            if name == "hex_terrain":
                cmap = ListedColormap(render.TERRAIN_COLORS)
                norm = BoundaryNorm(
                    np.arange(-0.5, len(render.TERRAIN_COLORS) + 0.5),
                    cmap.N,
                )
                ax.imshow(compiled.terrain, cmap=cmap, norm=norm, interpolation="nearest")
            else:
                raster = render.to_raster(grid, fields[name], width=180, height=90)
                if name == "elevation":
                    ax.imshow(raster, cmap=render.ELEVATION_CMAP,
                              norm=render.ELEVATION_NORM,
                              extent=[-180, 180, -90, 90])
                elif name == "crust_age":
                    ax.imshow(raster, cmap="cividis",
                              vmin=0.0,
                              vmax=max(float(np.nanpercentile(fields[name], 98)), 1.0),
                              extent=[-180, 180, -90, 90])
                else:
                    ax.imshow(raster, cmap="tab20", extent=[-180, 180, -90, 90])
            if row_i == 0:
                ax.set_title(name, fontsize=9)
            if col_i == 0:
                ax.set_ylabel(f"{int(row['cells'])} cells", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def promotion_audit_summary(
    *,
    stage_matrix: dict[str, Any],
    highres_assets: dict[str, Any],
    ci_review: dict[str, Any],
    ci_rows: list[dict[str, Any]],
    root_p90_inventory: dict[str, Any],
) -> dict[str, Any]:
    root_metrics = root_p90_inventory["metrics"]
    owner_counts = dict(root_p90_inventory.get("owner_counts", {}))
    non_asset_gap_count = int(sum(
        1
        for gap in root_p90_inventory.get("gaps", ())
        if gap.get("category") != "asset_review_pending"
    ))
    owner_blockers = tuple(
        f"{owner}_residuals_unresolved"
        for owner, count in sorted(owner_counts.items())
        if int(count) > 0 and owner != "compiler_render"
    )
    blockers: list[str] = []
    if stage_matrix["pass_count"] != stage_matrix["suite_count"]:
        blockers.append("p76_p90_stage_diagnostic_gate_failure")
    if not highres_assets["all_required_assets_present"]:
        blockers.append("highres_review_assets_incomplete")
    if not highres_assets.get("p69_highres_ready_all_members", False):
        blockers.append("p69_highres_members_not_ready")
    if not highres_assets.get("p69_promotion_ready", False):
        for blocker in highres_assets.get("p69_promotion_blockers", ()):
            blockers.append(f"p69_{blocker}")
    if not ci_review["all_required_assets_present"]:
        blockers.append("ci_review_assets_incomplete")
    if not ci_review["all_compilers_passed"]:
        blockers.append("ci_compiler_consistency_failure")
    if not ci_review["all_gap_inventories_ready"]:
        blockers.append("ci_gap_inventory_incomplete")
    if non_asset_gap_count > 0:
        blockers.append("p90_current_world_residuals_unresolved")
        blockers.extend(owner_blockers)
    blockers = list(dict.fromkeys(blockers))

    audit_completed = bool(
        stage_matrix["suite_count"] == len(STAGE_SUITES)
        and stage_matrix["missing_count"] == 0
        and highres_assets["available"]
        and highres_assets["all_required_assets_present"]
        and ci_review["world_count"] == 2
        and ci_review["all_required_assets_present"]
        and ci_review["contact_sheet_present"]
    )
    promotion_ready = bool(
        audit_completed
        and stage_matrix["pass_count"] == stage_matrix["suite_count"]
        and highres_assets.get("p69_promotion_ready", False)
        and ci_review["all_compilers_passed"]
        and ci_review["all_gap_inventories_ready"]
        and non_asset_gap_count == 0
        and not blockers
    )
    decision_recorded = bool((promotion_ready and not blockers) or blockers)

    return {
        "schema": SCHEMA,
        "status": (
            "integrated_real_earth_morphology_promotion_audit_ready"
            if audit_completed and decision_recorded
            else "integrated_real_earth_morphology_promotion_audit_incomplete"
        ),
        "promotion_hierarchy": PROMOTION_HIERARCHY,
        "stage_matrix": stage_matrix,
        "highres_review_assets": highres_assets,
        "ci_review": ci_review,
        "ci_worlds": tuple(
            {
                key: value
                for key, value in row.items()
                if key not in {"world", "compiled", "inventory"}
            }
            for row in ci_rows
        ),
        "root_p90_inventory_digest": str(root_p90_inventory.get("summary_digest", "")),
        "root_p90_owner_counts": owner_counts,
        "root_p90_category_counts": dict(
            root_p90_inventory.get("category_counts", {})),
        "root_p90_current_residual_items": tuple(
            root_p90_inventory.get("current_residual_items", ())),
        "promotion_decision": {
            "audit_completed": audit_completed,
            "promotion_ready": promotion_ready,
            "promotion_decision_recorded": decision_recorded,
            "keeps_default_off_until_named_residuals_resolved": bool(
                decision_recorded and not promotion_ready),
            "promotion_blockers": tuple(blockers),
            "next_recommended_entry": (
                "P92.production_residual_owner_repair_plan"
                if blockers else "default_earthlike_promotion"
            ),
        },
        "metrics": {
            "stage_suite_count": int(stage_matrix["suite_count"]),
            "stage_suite_pass_count": int(stage_matrix["pass_count"]),
            "stage_missing_count": int(stage_matrix["missing_count"]),
            "stage_failing_count": int(stage_matrix["failing_count"]),
            "highres_member_count": int(highres_assets["member_count"]),
            "highres_required_asset_files": int(
                highres_assets["expected_required_asset_files"]),
            "highres_actual_required_asset_files": int(
                highres_assets["actual_required_asset_files"]),
            "highres_contact_sheet_present_count": int(
                highres_assets["contact_sheet_present_count"]),
            "ci_world_count": int(ci_review["world_count"]),
            "ci_required_asset_files": int(ci_review["expected_required_asset_files"]),
            "ci_actual_required_asset_files": int(ci_review["actual_required_asset_files"]),
            "ci_asset_set_complete_count": int(ci_review["asset_set_complete_count"]),
            "ci_compiler_passed_count": int(ci_review["compiler_passed_count"]),
            "ci_inventory_ready_count": int(ci_review["inventory_ready_count"]),
            "ci_contact_sheet_present": bool(ci_review["contact_sheet_present"]),
            "p90_gap_count": int(root_metrics["gap_count"]),
            "p90_non_asset_gap_count": int(non_asset_gap_count),
            "p90_owner_layer_count": int(root_metrics["owner_layer_count"]),
            "p90_category_count": int(root_metrics["category_count"]),
            "p90_unassigned_gap_count": int(root_metrics["unassigned_gap_count"]),
            "p90_generic_blocker_count": int(root_metrics["generic_blocker_count"]),
            "p90_compiler_passed_envelope": bool(
                root_metrics["compiler_passed_envelope"]),
            "p90_high_flat_interior_fraction_of_continental_land": float(
                root_metrics["high_flat_interior_fraction_of_continental_land"]),
            "p90_basin_lowland_fraction_of_continental_land": float(
                root_metrics["basin_lowland_fraction_of_continental_land"]),
            "promotion_blocker_count": int(len(blockers)),
            "promotion_ready": bool(promotion_ready),
            "audit_completed": bool(audit_completed),
            "promotion_decision_recorded": bool(decision_recorded),
        },
        "acceptance": {
            "p76_p90_stage_summaries_available": bool(
                stage_matrix["missing_count"] == 0),
            "p76_p90_stage_gates_pass": bool(
                stage_matrix["pass_count"] == stage_matrix["suite_count"]),
            "p69_highres_review_assets_available": bool(
                highres_assets["all_required_assets_present"]),
            "p69_highres_contact_sheet_available": bool(
                highres_assets["contact_sheet_present_count"] >= 1),
            "ci_900_and_2500_worlds_audited": bool(
                set(ci_review["cells"]) == {900, 2500}),
            "ci_required_assets_available": bool(
                ci_review["all_required_assets_present"]),
            "ci_contact_sheet_available": bool(ci_review["contact_sheet_present"]),
            "ci_compiler_consistency_passed": bool(
                ci_review["all_compilers_passed"]),
            "ci_gap_inventories_ready": bool(
                ci_review["all_gap_inventories_ready"]),
            "p90_residuals_mapped": bool(
                root_metrics["gap_count"] > 0
                and root_metrics["unassigned_gap_count"] == 0
                and root_metrics["generic_blocker_count"] == 0),
            "promotion_decision_recorded": bool(decision_recorded),
            "keeps_default_off_until_named_residuals_resolved": bool(
                decision_recorded and not promotion_ready),
            "next_recommended_entry_defined": bool(
                blockers or promotion_ready),
        },
    }
