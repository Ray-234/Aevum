"""Regression gate for P26 continent/margin geometry work.

P26 changes are expensive to validate because the failure mode is world-level:
an isolated continent-widening fixture can pass while the full Earth-like run
regresses land fraction, ribbon geometry, ocean-basin partitioning, or seam
continuity.  This module compares existing P12 release summaries and records a
small, explicit no-regression contract against the accepted P25b baseline.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


SCHEMA = "aevum.p26_regression_gate.v1"


def _entry(summary: dict[str, Any], preset: str) -> dict[str, Any]:
    entries = summary.get("entries", [])
    for entry in entries:
        if str(entry.get("preset", "")) == preset:
            return entry
    raise KeyError(f"preset {preset!r} not present in P12 summary")


def _metric(entry: dict[str, Any], path: tuple[str, ...], default: float = 0.0) -> float:
    cur: Any = entry
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return float(default)
        cur = cur[key]
    try:
        return float(cur)
    except (TypeError, ValueError):
        return float(default)


def _validation_failures(entry: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    for name, check in entry.get("validation", {}).items():
        if not isinstance(check, dict):
            continue
        if check.get("passed", True):
            continue
        details = check.get("hard_failures", [])
        if details:
            failures.extend(f"{name}: {item}" for item in details)
        else:
            failures.append(str(name))
    failures.extend(str(item) for item in entry.get("release_gate", {}).get("failures", []))
    return failures


def _check_at_least(name: str, baseline: float, candidate: float, *,
                    min_value: float | None = None,
                    tolerance: float = 0.0,
                    meaning: str) -> dict[str, Any]:
    floor = baseline - tolerance
    if min_value is not None:
        floor = max(floor, min_value)
    return {
        "name": name,
        "passed": bool(candidate >= floor),
        "baseline": baseline,
        "candidate": candidate,
        "limit": floor,
        "direction": "at_least",
        "failure_meaning": meaning,
    }


def _check_at_most(name: str, baseline: float, candidate: float, *,
                   max_value: float | None = None,
                   tolerance: float = 0.0,
                   meaning: str) -> dict[str, Any]:
    ceiling = baseline + tolerance
    if max_value is not None:
        ceiling = min(ceiling, max_value)
    return {
        "name": name,
        "passed": bool(candidate <= ceiling),
        "baseline": baseline,
        "candidate": candidate,
        "limit": ceiling,
        "direction": "at_most",
        "failure_meaning": meaning,
    }


def _check_range(name: str, baseline: float, candidate: float, *,
                 lower_ratio: float,
                 upper_ratio: float,
                 absolute_floor: float,
                 meaning: str) -> dict[str, Any]:
    lower = max(absolute_floor, baseline * lower_ratio)
    upper = max(lower, baseline * upper_ratio)
    return {
        "name": name,
        "passed": bool(lower <= candidate <= upper),
        "baseline": baseline,
        "candidate": candidate,
        "limit": [lower, upper],
        "direction": "range",
        "failure_meaning": meaning,
    }


def _check_ocean_basin_count_no_regression(
    name: str,
    baseline: float,
    candidate: float,
    *,
    meaning: str,
) -> dict[str, Any]:
    """Avoid locking in under-partitioned ocean basins as the regression target."""
    ratio_lower = max(2.0, baseline * 0.55)
    ratio_upper = max(ratio_lower, baseline * 1.75)
    earthlike_lower = 3.0
    earthlike_upper = 16.0
    ratio_pass = ratio_lower <= candidate <= ratio_upper
    earthlike_pass = earthlike_lower <= candidate <= earthlike_upper
    return {
        "name": name,
        "passed": bool(ratio_pass or earthlike_pass),
        "baseline": baseline,
        "candidate": candidate,
        "limit": {
            "ratio_range": [ratio_lower, ratio_upper],
            "earthlike_range": [earthlike_lower, earthlike_upper],
        },
        "direction": "range",
        "failure_meaning": meaning,
    }


def compare_p12_summaries(
    baseline: dict[str, Any],
    candidate: dict[str, Any],
    *,
    preset: str = "earthlike",
) -> dict[str, Any]:
    """Compare a candidate P12 summary against an accepted baseline summary."""
    base_entry = _entry(baseline, preset)
    cand_entry = _entry(candidate, preset)
    checks: list[dict[str, Any]] = []

    failures = _validation_failures(cand_entry)
    checks.append({
        "name": "release_and_validation_pass",
        "passed": not failures,
        "baseline": baseline.get("release_decision", {}).get("status", ""),
        "candidate": candidate.get("release_decision", {}).get("status", ""),
        "limit": "no release failures and no validation hard failures",
        "direction": "boolean",
        "failure_meaning": "; ".join(failures[:6]) if failures else "",
    })

    land_base = _metric(base_entry, ("land_fraction",))
    land_cand = _metric(cand_entry, ("land_fraction",))
    checks.append(_check_at_least(
        "land_fraction_no_regression",
        land_base,
        land_cand,
        min_value=0.25 if land_base >= 0.25 else None,
        tolerance=0.003,
        meaning="P26 must not undo the P25b land-budget repair.",
    ))

    comp_base = _metric(base_entry, ("morphology", "land_component_count"))
    comp_cand = _metric(cand_entry, ("morphology", "land_component_count"))
    checks.append(_check_at_most(
        "land_component_count_no_regression",
        comp_base,
        comp_cand,
        max_value=14.0,
        tolerance=2.0,
        meaning="P26 must not re-fragment exposed land into many islands/components.",
    ))

    ribbon_base = _metric(base_entry, ("morphology", "land_ribbon_fraction_gt_0_5"))
    ribbon_cand = _metric(cand_entry, ("morphology", "land_ribbon_fraction_gt_0_5"))
    checks.append(_check_at_most(
        "land_ribbon_no_regression",
        ribbon_base,
        ribbon_cand,
        max_value=0.55,
        tolerance=0.015,
        meaning="P26 must not worsen the main Earth-reference ribbon problem.",
    ))

    cont_ribbon_base = _metric(
        base_entry, ("morphology", "continental_ribbon_fraction_gt_0_5"))
    cont_ribbon_cand = _metric(
        cand_entry, ("morphology", "continental_ribbon_fraction_gt_0_5"))
    checks.append(_check_at_most(
        "continental_ribbon_no_regression",
        cont_ribbon_base,
        cont_ribbon_cand,
        tolerance=0.015,
        meaning="P26 upstream geometry must not make continental crust narrower.",
    ))

    coast_base = _metric(base_entry, ("morphology", "land_coastline_complexity_largest"))
    coast_cand = _metric(cand_entry, ("morphology", "land_coastline_complexity_largest"))
    checks.append(_check_at_most(
        "largest_coastline_complexity_no_regression",
        coast_base,
        coast_cand,
        tolerance=1.0,
        meaning="P26 must not increase largest-landmass coastline jaggedness.",
    ))

    basin_base = _metric(base_entry, ("ocean_geography", "basin_count"))
    basin_cand = _metric(cand_entry, ("ocean_geography", "basin_count"))
    checks.append(_check_ocean_basin_count_no_regression(
        "ocean_basin_count_no_regression",
        basin_base,
        basin_cand,
        meaning="P26 must not collapse or over-fragment ocean-basin partitioning.",
    ))

    seam_basin = _metric(cand_entry, ("seam_continuity", "seam_ocean_basin_mismatch_fraction"))
    seam_render = _metric(
        cand_entry, ("seam_continuity", "render_duplicate_basin_mismatch_fraction"))
    checks.append({
        "name": "ocean_basin_seam_continuity",
        "passed": bool(seam_basin == 0.0 and seam_render == 0.0),
        "baseline": 0.0,
        "candidate": {
            "seam_ocean_basin_mismatch_fraction": seam_basin,
            "render_duplicate_basin_mismatch_fraction": seam_render,
        },
        "limit": 0.0,
        "direction": "equals",
        "failure_meaning": "P26 must preserve spherical/dateline ocean-basin continuity.",
    })

    failed = [check for check in checks if not check["passed"]]
    return {
        "schema": SCHEMA,
        "preset": preset,
        "status": "pass" if not failed else "fail",
        "checks": checks,
        "failed_checks": [check["name"] for check in failed],
        "acceptance": {
            "no_world_level_regressions": not failed,
            "release_and_validation_pass": checks[0]["passed"],
            "land_planform_not_regressed": all(
                check["passed"]
                for check in checks
                if check["name"] in {
                    "land_fraction_no_regression",
                    "land_component_count_no_regression",
                    "land_ribbon_no_regression",
                    "continental_ribbon_no_regression",
                    "largest_coastline_complexity_no_regression",
                }
            ),
            "ocean_partition_not_regressed": all(
                check["passed"]
                for check in checks
                if check["name"] in {
                    "ocean_basin_count_no_regression",
                    "ocean_basin_seam_continuity",
                }
            ),
        },
    }


def compare_summary_files(
    baseline_path: Path,
    candidate_path: Path,
    *,
    preset: str = "earthlike",
) -> dict[str, Any]:
    baseline = json.loads(Path(baseline_path).read_text())
    candidate = json.loads(Path(candidate_path).read_text())
    report = compare_p12_summaries(baseline, candidate, preset=preset)
    report["baseline_path"] = str(baseline_path)
    report["candidate_path"] = str(candidate_path)
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m aevum.diagnostics.p26_regression_gate",
        description="Compare candidate P12 summary against a P25b/P26 baseline.",
    )
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--preset", default="earthlike")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    report = compare_summary_files(args.baseline, args.candidate, preset=args.preset)
    text = json.dumps(report, indent=2) + "\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print("== aevum :: P26 regression gate ==")
    print(f"   status: {report['status']}")
    print(f"   failed checks: {', '.join(report['failed_checks']) if report['failed_checks'] else 'none'}")
    if args.out is not None:
        print(f"   wrote {args.out}")
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
