"""Biome-envelope gate against Koppen proxy and RESOLVE references."""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "aevum.earth_climate_biome_gate.v1"
BIOME_LABELS = {
    0: "ocean",
    1: "ice",
    2: "desert",
    3: "grassland",
    4: "forest",
    5: "tundra",
    6: "tropical",
}
RESOLVE_TO_AEVUM = {
    0: 0,   # no-data / ocean
    1: 6,   # tropical moist broadleaf
    2: 6,   # tropical dry broadleaf
    3: 6,   # tropical coniferous
    4: 4,   # temperate broadleaf/mixed
    5: 4,   # temperate conifer
    6: 4,   # boreal forest / taiga
    7: 3,   # tropical grassland/savanna
    8: 3,   # temperate grassland
    9: 3,   # flooded grassland
    10: 3,  # montane grassland
    11: 5,  # tundra
    12: 4,  # Mediterranean woodland/scrub
    13: 2,  # desert/xeric shrubland
    14: 6,  # mangrove
}


@dataclass(frozen=True)
class EarthClimateBiomeGateConfig:
    earth_reference_npz: Path
    terminal_summary_json: Path
    outdir: Path


def _json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _area_fraction(area: np.ndarray, mask: np.ndarray, denom_mask: np.ndarray) -> float:
    area = np.asarray(area, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    denom_mask = np.asarray(denom_mask, dtype=bool)
    denom = float(np.nansum(area[denom_mask]))
    if denom <= 0.0:
        return float("nan")
    return float(np.nansum(area[mask & denom_mask]) / denom)


def _map_resolve_to_aevum(resolve: np.ndarray) -> np.ndarray:
    resolve = np.asarray(resolve, dtype=np.int16)
    out = np.zeros(resolve.shape, dtype=np.int16)
    for src, dst in RESOLVE_TO_AEVUM.items():
        out[resolve == src] = dst
    return out


def _biome_metrics(
    label: str,
    preset: str,
    seed: int,
    area: np.ndarray,
    land: np.ndarray,
    biome: np.ndarray,
) -> dict[str, Any]:
    biome = np.asarray(biome, dtype=np.int16)
    land = np.asarray(land, dtype=bool)
    row: dict[str, Any] = {
        "label": label,
        "preset": preset,
        "seed": seed,
        "land_fraction": _area_fraction(area, land, np.ones_like(land, dtype=bool)),
        "forest_tropical_land_fraction": _area_fraction(
            area, land & ((biome == 4) | (biome == 6)), land),
        "tropical_land_fraction": _area_fraction(area, land & (biome == 6), land),
        "forest_land_fraction": _area_fraction(area, land & (biome == 4), land),
        "desert_land_fraction": _area_fraction(area, land & (biome == 2), land),
        "grassland_land_fraction": _area_fraction(area, land & (biome == 3), land),
        "tundra_ice_land_fraction": _area_fraction(
            area, land & ((biome == 1) | (biome == 5)), land),
    }
    for code, name in BIOME_LABELS.items():
        row[f"{name}_land_fraction"] = _area_fraction(area, land & (biome == code), land)
    return row


def _earth_metrics(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as z:
        area = np.asarray(z["cell_area"], dtype=np.float64)
        land = np.asarray(z["earth__land_mask"], dtype=bool)
        proxy = np.asarray(z["earth__biome_class_proxy"], dtype=np.int16)
        resolve = np.asarray(z["earth__resolve_biome_class"], dtype=np.int16)
    resolve_coarse = _map_resolve_to_aevum(resolve)
    return {
        "koppen_proxy": _biome_metrics(
            "Earth Koppen proxy", "earth_reference", -1, area, land, proxy),
        "resolve_coarse": _biome_metrics(
            "Earth RESOLVE coarse", "earth_reference", -1,
            area, land & (resolve > 0), resolve_coarse),
    }


def _generated_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    arrays = Path(summary["arrays"])
    with np.load(arrays, allow_pickle=False) as z:
        area = np.asarray(z["cell_area"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        land = np.asarray(z["terrain__elevation_m"], dtype=np.float64) >= sea
        biome = np.asarray(z["biosphere__biome"], dtype=np.float64).astype(np.int16)
    row = _biome_metrics(
        Path(str(summary.get("assets_dir", ""))).name,
        str(summary.get("preset", "")),
        int(summary.get("seed", -1)),
        area,
        land,
        biome,
    )
    row["arrays"] = str(arrays)
    row["mode"] = (
        "earthlike_calibration"
        if "earthlike" in row["preset"].lower()
        else "diagnostic_only"
    )
    return row


def _ratio(value: float, ref: float) -> float:
    if not np.isfinite(value) or not np.isfinite(ref) or abs(ref) <= 1.0e-12:
        return float("nan")
    return float(value / ref)


def _check(
    label: str,
    group: str,
    metric: str,
    generated: float,
    reference: float,
    operator: str,
    threshold: float,
    severity: str,
    message: str,
) -> dict[str, Any]:
    ratio = _ratio(generated, reference)
    skipped = not np.isfinite(float(generated))
    if skipped:
        passed = True
    elif operator == "ratio>=":
        passed = np.isfinite(ratio) and ratio >= threshold
    elif operator == "ratio<=":
        passed = np.isfinite(ratio) and ratio <= threshold
    elif operator == ">=":
        passed = float(generated) >= threshold
    elif operator == "<=":
        passed = float(generated) <= threshold
    else:
        raise ValueError(f"unsupported operator {operator!r}")
    return {
        "label": label,
        "group": group,
        "metric": metric,
        "generated": float(generated),
        "reference": float(reference),
        "ratio_to_reference": ratio,
        "operator": operator,
        "threshold": float(threshold),
        "severity": severity,
        "passed": bool(passed),
        "skipped": bool(skipped),
        "message": message,
    }


def _earthlike_checks(
    row: dict[str, Any],
    references: dict[str, Any],
) -> list[dict[str, Any]]:
    proxy = references["koppen_proxy"]
    resolve = references["resolve_coarse"]
    label = str(row["label"])
    return [
        _check(
            label,
            "koppen_proxy",
            "forest_tropical_land_fraction",
            row["forest_tropical_land_fraction"],
            proxy["forest_tropical_land_fraction"],
            "ratio>=",
            0.40,
            "fail",
            "forest+tropical biome envelope is low versus Koppen-derived proxy",
        ),
        _check(
            label,
            "koppen_proxy",
            "tropical_land_fraction",
            row["tropical_land_fraction"],
            proxy["tropical_land_fraction"],
            "ratio>=",
            0.35,
            "fail",
            "tropical biome envelope is low versus Koppen-derived proxy",
        ),
        _check(
            label,
            "resolve_coarse",
            "forest_tropical_land_fraction",
            row["forest_tropical_land_fraction"],
            resolve["forest_tropical_land_fraction"],
            "ratio>=",
            0.40,
            "warn",
            "forest+tropical biome envelope is low versus RESOLVE coarse grouping",
        ),
        _check(
            label,
            "koppen_proxy",
            "desert_land_fraction",
            row["desert_land_fraction"],
            proxy["desert_land_fraction"],
            "ratio<=",
            2.25,
            "warn",
            "desert biome envelope should not greatly exceed Earth proxy",
        ),
        _check(
            label,
            "koppen_proxy",
            "tundra_ice_land_fraction",
            row["tundra_ice_land_fraction"],
            proxy["tundra_ice_land_fraction"],
            "ratio>=",
            0.50,
            "warn",
            "cold biome envelope is low versus Earth proxy",
        ),
    ]


def _write_csv(path: Path, rows: list[dict[str, Any]], keys: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in keys})
    return path


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Biome Gate",
        "",
        f"Verdict: `{report['verdict']}`",
        f"Failures: `{report['failure_count']}`",
        f"Warnings: `{report['warning_count']}`",
        "",
        "## Reference Envelopes",
        "",
    ]
    for name, ref in report["earth_references"].items():
        lines.extend([
            f"### {name}",
            "",
            f"- Forest+tropical land fraction: `{ref['forest_tropical_land_fraction']:.3f}`",
            f"- Tropical land fraction: `{ref['tropical_land_fraction']:.3f}`",
            f"- Desert land fraction: `{ref['desert_land_fraction']:.3f}`",
            f"- Tundra/ice land fraction: `{ref['tundra_ice_land_fraction']:.3f}`",
            "",
        ])
    lines.extend(["## Failed / Warning Checks", ""])
    for row in report["checks"]:
        if row["passed"]:
            continue
        status = "FAIL" if row["severity"] == "fail" else "WARN"
        lines.append(
            f"- {status} `{row['label']}` `{row['group']}` `{row['metric']}` "
            f"generated `{row['generated']:.3f}`, ref `{row['reference']:.3f}`, "
            f"ratio `{row['ratio_to_reference']:.3f}`: {row['message']}"
        )
    lines.extend(["", "## Generated Earthlike Biomes", ""])
    for row in report["generated_metrics"]:
        if row.get("mode") != "earthlike_calibration":
            continue
        lines.extend([
            f"### {row['label']}",
            "",
            f"- Forest+tropical: `{row['forest_tropical_land_fraction']:.3f}`",
            f"- Forest / tropical: `{row['forest_land_fraction']:.3f}` / `{row['tropical_land_fraction']:.3f}`",
            f"- Desert / grassland: `{row['desert_land_fraction']:.3f}` / `{row['grassland_land_fraction']:.3f}`",
            f"- Tundra+ice: `{row['tundra_ice_land_fraction']:.3f}`",
            "",
        ])
    return "\n".join(lines)


def run_earth_climate_biome_gate(
    config: EarthClimateBiomeGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    terminal = json.loads(Path(config.terminal_summary_json).read_text())
    refs = _earth_metrics(Path(config.earth_reference_npz))
    generated = [_generated_metrics(row) for row in terminal.get("summaries", [])]
    generated.sort(key=lambda row: (
        row.get("mode") != "earthlike_calibration",
        str(row.get("label")),
    ))
    checks: list[dict[str, Any]] = []
    for row in generated:
        if row["mode"] == "earthlike_calibration":
            checks.extend(_earthlike_checks(row, refs))
    failures = [row for row in checks if not row["passed"] and row["severity"] == "fail"]
    warnings = [row for row in checks if not row["passed"] and row["severity"] == "warn"]
    skipped = [row for row in checks if row.get("skipped")]
    if failures:
        verdict = "fail"
    elif warnings:
        verdict = "pass_with_warnings"
    else:
        verdict = "pass"

    metric_keys = [
        "label", "preset", "seed", "mode", "land_fraction",
        "forest_tropical_land_fraction", "forest_land_fraction",
        "tropical_land_fraction", "desert_land_fraction",
        "grassland_land_fraction", "tundra_ice_land_fraction",
    ]
    reference_rows = [
        refs["koppen_proxy"] | {"mode": "earth_reference_koppen_proxy"},
        refs["resolve_coarse"] | {"mode": "earth_reference_resolve_coarse"},
    ]
    metrics_csv = _write_csv(
        outdir / "earth_climate_biome_metrics.csv",
        reference_rows + generated,
        metric_keys,
    )
    checks_csv = _write_csv(
        outdir / "earth_climate_biome_checks.csv",
        checks,
        [
            "label", "group", "metric", "generated", "reference",
            "ratio_to_reference", "operator", "threshold", "severity",
            "passed", "skipped", "message",
        ],
    )
    report = {
        "schema": SCHEMA,
        "earth_reference_npz": str(config.earth_reference_npz),
        "terminal_summary_json": str(config.terminal_summary_json),
        "verdict": verdict,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "earth_references": refs,
        "generated_metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "skipped": skipped,
        "metrics_csv": str(metrics_csv),
        "checks_csv": str(checks_csv),
    }
    md_path = outdir / "earth_climate_biome_gate_report.md"
    report["markdown"] = str(md_path)
    (outdir / "earth_climate_biome_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default)
    )
    md_path.write_text(_render_markdown(report))
    return report
