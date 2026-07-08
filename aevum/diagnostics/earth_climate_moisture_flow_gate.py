"""C4e seasonal moisture-flow-network gate.

This gate evaluates the diagnostic object layer added after C4d regional
hydroclimate.  It intentionally does not score total precipitation: scalar and
Earth-pattern climate gates still own that contract.  C4e is judged on whether
seasonal source-ocean moisture, downwind land pathways, and flow-network
objects are archived and read as coherent routed corridors.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

from aevum.diagnostics.earth_climate_hydro_region_gate import (
    _array_path,
    _check,
    _component_shape_metrics,
    _corr,
    _edges_from_latlon,
    _json_default,
    _prefixed,
    _safe_float,
    _seasonal_map_metrics,
    _values,
    _write_csv,
)


SCHEMA = "aevum.earth_climate_moisture_flow_gate.v1"


@dataclass(frozen=True)
class EarthClimateMoistureFlowGateConfig:
    terminal_summary_json: Path
    outdir: Path
    render_contact_sheets: bool = True


def _candidate_network_paths(row: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    raw = row.get("moisture_flow_networks_json")
    if raw:
        paths.append(Path(str(raw)))
    assets = row.get("assets_dir")
    if assets:
        paths.append(Path(str(assets)) / "moisture_flow_networks.json")
    return paths


def _load_networks(row: dict[str, Any]) -> tuple[list[dict[str, Any]], Path | None]:
    for path in _candidate_network_paths(row):
        if path.exists():
            payload = json.loads(path.read_text())
            networks = payload.get("networks", [])
            if isinstance(networks, list):
                return [obj for obj in networks if isinstance(obj, dict)], path
            return [], path
    return [], None


def _max(rows: list[dict[str, Any]], key: str) -> float:
    values = _values(rows, key)
    return float(max(values)) if values else 0.0


def _sum(rows: list[dict[str, Any]], key: str) -> float:
    values = _values(rows, key)
    return float(sum(values)) if values else 0.0


def _pct(rows: list[dict[str, Any]], key: str, q: float) -> float:
    values = _values(rows, key)
    return float(np.percentile(values, q)) if values else 0.0


def _kind(rows: list[dict[str, Any]], kind: str) -> list[dict[str, Any]]:
    return [row for row in rows if str(row.get("kind", "")) == kind]


def _season_count(rows: list[dict[str, Any]]) -> int:
    return int(len({str(row.get("season", "")) for row in rows if row.get("season")}))


def _median(values: list[float]) -> float:
    finite = [float(value) for value in values if np.isfinite(value)]
    return float(np.median(finite)) if finite else 0.0


def _seasonal_binary_map_metrics(
    active4: np.ndarray,
    domain: np.ndarray,
    area: np.ndarray,
    edges: np.ndarray,
) -> dict[str, float]:
    active4 = np.asarray(active4, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    area = np.asarray(area, dtype=np.float64)
    if active4.ndim != 2 or active4.shape[0] != 4 or active4.shape[1] != domain.size:
        return {
            "component_count_p50": 0.0,
            "largest_component_share_p50": 0.0,
            "active_world_fraction_p50": 0.0,
            "active_domain_fraction_p50": 0.0,
            "boundary_per_active_cell_p50": 0.0,
        }
    domain_area = max(float(np.sum(area[domain])), 1.0e-12)
    rows: list[dict[str, float]] = []
    for season in range(4):
        active = active4[season] & domain
        metrics = _component_shape_metrics(active, area, edges)
        metrics["active_domain_fraction"] = float(np.sum(area[active]) / domain_area)
        rows.append(metrics)

    def med(key: str) -> float:
        return _median([row[key] for row in rows])

    return {
        "component_count_p50": med("component_count"),
        "largest_component_share_p50": med("largest_component_share"),
        "active_world_fraction_p50": med("active_world_fraction"),
        "active_domain_fraction_p50": med("active_domain_fraction"),
        "boundary_per_active_cell_p50": med("boundary_per_active_cell"),
    }


def _seasonal_domain_pctl(field: np.ndarray, domain: np.ndarray, q: float) -> float:
    field = np.asarray(field, dtype=np.float64)
    domain = np.asarray(domain, dtype=bool)
    if field.ndim != 2 or field.shape[0] != 4 or field.shape[1] != domain.size:
        return 0.0
    values: list[float] = []
    for season in range(4):
        vals = field[season, domain]
        vals = vals[np.isfinite(vals)]
        values.append(float(np.percentile(vals, q)) if vals.size else 0.0)
    return _median(values)


def _seasonal_corr_median(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    if a.ndim != 2 or b.shape != a.shape or a.shape[0] != 4:
        return 0.0
    return _median([_corr(a[season], b[season], mask) for season in range(4)])


def _array_metrics(summary_row: dict[str, Any]) -> dict[str, Any]:
    defaults = {
        "arrays_found": 0.0,
        "arrays_path": "",
        "flow_map_edges_found": 0.0,
        "moisture_source_ocean_p90_median": 0.0,
        "moisture_pathway_land_p90_median": 0.0,
        "moisture_pathway_ocean_p90_median": 0.0,
        "pathway_support_corr_median": 0.0,
        "pathway_precip_corr_median": 0.0,
        "pathway_rain_shadow_corr_median": 0.0,
        "pathway_map_component_count_p50": 0.0,
        "pathway_map_largest_component_share_p50": 0.0,
        "pathway_map_active_world_fraction_p50": 0.0,
        "pathway_map_active_land_fraction_p50": 0.0,
        "pathway_map_boundary_per_active_cell_p50": 0.0,
        "network_id_map_component_count_p50": 0.0,
        "network_id_map_largest_component_share_p50": 0.0,
        "network_id_map_active_world_fraction_p50": 0.0,
        "network_id_map_active_domain_fraction_p50": 0.0,
        "network_id_map_boundary_per_active_cell_p50": 0.0,
    }
    path = _array_path(summary_row)
    if path is None:
        return defaults

    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "cell_area",
            "terrain__elevation_m",
            "sea_level_m",
            "atmosphere__moisture_flow_source",
            "atmosphere__moisture_flow_pathway",
            "climate__moisture_flow_network_id",
            "climate__seasonal_precipitation",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
            "climate__rain_shadow_index",
        )
        if any(key not in z for key in required):
            out = dict(defaults)
            out["arrays_found"] = 1.0
            out["arrays_path"] = str(path)
            return out
        lat = np.asarray(z["lat"], dtype=np.float64)
        lon = np.asarray(z["lon"], dtype=np.float64)
        area = np.asarray(z["cell_area"], dtype=np.float64)
        elev = np.asarray(z["terrain__elevation_m"], dtype=np.float64)
        sea = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        source = np.asarray(z["atmosphere__moisture_flow_source"], dtype=np.float64)
        pathway = np.asarray(z["atmosphere__moisture_flow_pathway"], dtype=np.float64)
        network_id = np.asarray(z["climate__moisture_flow_network_id"], dtype=np.float64)
        precip = np.asarray(z["climate__seasonal_precipitation"], dtype=np.float64)
        monsoon = np.asarray(z["climate__monsoon_rainfall_corridor"], dtype=np.float64)
        storm = np.asarray(z["climate__storm_track_rainfall_corridor"], dtype=np.float64)
        shadow = np.asarray(z["climate__rain_shadow_index"], dtype=np.float64)

    land = elev >= sea
    ocean = ~land
    edges = _edges_from_latlon(lat, lon)
    support = np.maximum(monsoon, storm)
    network_active = np.isfinite(network_id) & (network_id > 0.0)

    out = {
        "arrays_found": 1.0,
        "arrays_path": str(path),
        "flow_map_edges_found": 1.0 if edges.size else 0.0,
        "moisture_source_ocean_p90_median": _seasonal_domain_pctl(source, ocean, 90.0),
        "moisture_pathway_land_p90_median": _seasonal_domain_pctl(pathway, land, 90.0),
        "moisture_pathway_ocean_p90_median": _seasonal_domain_pctl(pathway, ocean, 90.0),
        "pathway_support_corr_median": _seasonal_corr_median(pathway, support, land),
        "pathway_precip_corr_median": _seasonal_corr_median(pathway, precip, land),
        "pathway_rain_shadow_corr_median": _seasonal_corr_median(pathway, shadow, land),
    }
    out.update(_prefixed(
        "pathway_map",
        _seasonal_map_metrics(pathway, land, area, edges, floor=0.035, percentile=82.0),
    ))
    out.update(_prefixed(
        "network_id_map",
        _seasonal_binary_map_metrics(network_active, land, area, edges),
    ))
    return out


def _generated_row(summary_row: dict[str, Any]) -> dict[str, Any]:
    label = Path(str(summary_row.get("assets_dir", summary_row.get("label", "")))).name
    networks, path = _load_networks(summary_row)
    kinds = Counter(str(row.get("kind", "unknown")) for row in networks)
    seasons = Counter(str(row.get("season", "unknown")) for row in networks)
    monsoon = _kind(networks, "monsoon_moisture_flow_network")
    storm = _kind(networks, "storm_track_moisture_flow_network")
    mixed = _kind(networks, "mixed_moisture_flow_network")
    cell_counts = _values(networks, "cell_count")
    row = {
        "label": label,
        "preset": str(summary_row.get("preset", "")),
        "seed": int(summary_row.get("seed", -1)),
        "assets_dir": str(summary_row.get("assets_dir", "")),
        "networks_json": str(path) if path else "",
        "object_json_found": 1.0 if path is not None else 0.0,
        "object_count": int(len(networks)),
        "kind_count": int(len(kinds)),
        "season_count": int(len(seasons)),
        "kind_counts": dict(sorted(kinds.items())),
        "season_counts": dict(sorted(seasons.items())),
        "largest_area_fraction": _max(networks, "area_fraction"),
        "top5_area_fraction_sum": float(sum(sorted(
            _values(networks, "area_fraction"), reverse=True)[:5])),
        "cell_count_p50": float(np.percentile(cell_counts, 50)) if cell_counts else 0.0,
        "cell_count_p90": float(np.percentile(cell_counts, 90)) if cell_counts else 0.0,
        "mean_pathway_p50": _pct(networks, "mean_pathway", 50),
        "mean_pathway_p90": _pct(networks, "mean_pathway", 90),
        "mean_precipitation_p50": _pct(networks, "mean_precipitation", 50),
        "monsoon_flow_object_count": int(len(monsoon)),
        "storm_track_flow_object_count": int(len(storm)),
        "mixed_flow_object_count": int(len(mixed)),
        "monsoon_flow_season_count": _season_count(monsoon),
        "storm_track_flow_season_count": _season_count(storm),
        "largest_monsoon_flow_area_fraction": _max(monsoon, "area_fraction"),
        "largest_storm_track_flow_area_fraction": _max(storm, "area_fraction"),
        "largest_mixed_flow_area_fraction": _max(mixed, "area_fraction"),
        "monsoon_flow_area_fraction_sum": _sum(monsoon, "area_fraction"),
        "storm_track_flow_area_fraction_sum": _sum(storm, "area_fraction"),
        "mixed_flow_area_fraction_sum": _sum(mixed, "area_fraction"),
    }
    row.update(_array_metrics(summary_row))
    return row


def _earthlike_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike runs must archive C4e moisture-flow objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="earthlike runs must archive C4e source/pathway arrays"),
        _check(label=label, group="map_archive", metric="flow_map_edges_found",
               generated=row["flow_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="C4e arrays must support map-readability checks"),
        _check(label=label, group="object_coverage", metric="object_count",
               generated=row["object_count"], operator=">=", threshold=25.0,
               severity="fail", message="earthlike worlds should expose many seasonal flow networks"),
        _check(label=label, group="object_coverage", metric="kind_count",
               generated=row["kind_count"], operator=">=", threshold=2.0,
               severity="fail", message="earthlike flow networks should include multiple moisture regimes"),
        _check(label=label, group="seasonal_coverage", metric="season_count",
               generated=row["season_count"], operator=">=", threshold=4.0,
               severity="fail", message="flow networks should be seasonal, not annual-only"),
        _check(label=label, group="monsoon_flow", metric="monsoon_flow_object_count",
               generated=row["monsoon_flow_object_count"], operator=">=", threshold=2.0,
               severity="fail", message="earthlike worlds should retain monsoon moisture-flow networks"),
        _check(label=label, group="storm_track_flow", metric="storm_track_flow_object_count",
               generated=row["storm_track_flow_object_count"], operator=">=", threshold=20.0,
               severity="fail", message="earthlike worlds should retain storm-track moisture-flow networks"),
        _check(label=label, group="coherence_proxy", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator=">=", threshold=0.010,
               severity="fail", message="flow networks should not be only tiny speckles"),
        _check(label=label, group="coherence_proxy", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator="<=", threshold=0.120,
               severity="warn", message="one flow network should not swallow a continental climate zone"),
        _check(label=label, group="source_pathway_strength",
               metric="moisture_source_ocean_p90_median",
               generated=row["moisture_source_ocean_p90_median"],
               operator=">=", threshold=0.25, severity="fail",
               message="source-ocean moisture should be visible"),
        _check(label=label, group="source_pathway_strength",
               metric="moisture_pathway_land_p90_median",
               generated=row["moisture_pathway_land_p90_median"],
               operator=">=", threshold=0.35, severity="fail",
               message="land moisture pathways should be visible"),
        _check(label=label, group="pathway_coupling",
               metric="pathway_support_corr_median",
               generated=row["pathway_support_corr_median"], operator=">=", threshold=0.10,
               severity="fail", message="flow pathways should align with monsoon/storm support"),
        _check(label=label, group="pathway_map_readability",
               metric="pathway_map_active_land_fraction_p50",
               generated=row["pathway_map_active_land_fraction_p50"],
               operator=">=", threshold=0.18, severity="fail",
               message="pathway maps should not be too sparse to read"),
        _check(label=label, group="pathway_map_readability",
               metric="pathway_map_active_land_fraction_p50",
               generated=row["pathway_map_active_land_fraction_p50"],
               operator="<=", threshold=0.45, severity="warn",
               message="pathway maps should remain corridor-like instead of wetting most land"),
        _check(label=label, group="pathway_map_readability",
               metric="pathway_map_largest_component_share_p50",
               generated=row["pathway_map_largest_component_share_p50"],
               operator=">=", threshold=0.14, severity="fail",
               message="pathways should include connected downwind corridors"),
        _check(label=label, group="pathway_map_readability",
               metric="pathway_map_boundary_per_active_cell_p50",
               generated=row["pathway_map_boundary_per_active_cell_p50"],
               operator="<=", threshold=2.80, severity="fail",
               message="pathways should not read as checkerboard texture"),
        _check(label=label, group="network_id_map_readability",
               metric="network_id_map_active_domain_fraction_p50",
               generated=row["network_id_map_active_domain_fraction_p50"],
               operator=">=", threshold=0.04, severity="fail",
               message="network-id maps should expose visible land corridors"),
        _check(label=label, group="network_id_map_readability",
               metric="network_id_map_largest_component_share_p50",
               generated=row["network_id_map_largest_component_share_p50"],
               operator=">=", threshold=0.10, severity="fail",
               message="network-id maps should not fragment into isolated cells"),
    ]


def _waterworld_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld diagnostics should archive C4e flow objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld diagnostics should archive C4e arrays"),
        _check(label=label, group="map_archive", metric="flow_map_edges_found",
               generated=row["flow_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="waterworld arrays must support C4e map checks"),
        _check(label=label, group="waterworld_false_positive", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator="<=", threshold=0.012,
               severity="fail", message="waterworld flow networks should remain island-scale"),
        _check(label=label, group="waterworld_false_positive",
               metric="largest_monsoon_flow_area_fraction",
               generated=row["largest_monsoon_flow_area_fraction"],
               operator="<=", threshold=0.008, severity="fail",
               message="waterworlds should not grow continent-scale monsoon flow networks"),
        _check(label=label, group="waterworld_false_positive",
               metric="pathway_map_active_world_fraction_p50",
               generated=row["pathway_map_active_world_fraction_p50"],
               operator="<=", threshold=0.035, severity="fail",
               message="waterworld pathway maps should stay small in world-area terms"),
        _check(label=label, group="waterworld_false_positive",
               metric="network_id_map_active_world_fraction_p50",
               generated=row["network_id_map_active_world_fraction_p50"],
               operator="<=", threshold=0.020, severity="fail",
               message="waterworld network-id maps should remain island-scale"),
    ]


def _arid_checks(row: dict[str, Any]) -> list[dict[str, Any]]:
    label = str(row["label"])
    return [
        _check(label=label, group="object_archive", metric="object_json_found",
               generated=row["object_json_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid diagnostics should archive C4e flow objects"),
        _check(label=label, group="array_archive", metric="arrays_found",
               generated=row["arrays_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid diagnostics should archive C4e arrays"),
        _check(label=label, group="map_archive", metric="flow_map_edges_found",
               generated=row["flow_map_edges_found"], operator=">=", threshold=1.0,
               severity="fail", message="arid arrays must support C4e map checks"),
        _check(label=label, group="object_coverage", metric="object_count",
               generated=row["object_count"], operator=">=", threshold=10.0,
               severity="fail", message="arid worlds should still expose routed moisture-flow diagnostics"),
        _check(label=label, group="arid_false_monsoon_guard",
               metric="largest_monsoon_flow_area_fraction",
               generated=row["largest_monsoon_flow_area_fraction"],
               operator="<=", threshold=0.040, severity="fail",
               message="arid worlds should not become broadly monsoonal"),
        _check(label=label, group="arid_flow_bounds", metric="largest_area_fraction",
               generated=row["largest_area_fraction"], operator="<=", threshold=0.180,
               severity="fail", message="arid flow networks should remain bounded"),
        _check(label=label, group="source_pathway_strength",
               metric="moisture_source_ocean_p90_median",
               generated=row["moisture_source_ocean_p90_median"],
               operator=">=", threshold=0.20, severity="fail",
               message="arid worlds should still have source-ocean moisture diagnostics"),
        _check(label=label, group="source_pathway_strength",
               metric="moisture_pathway_land_p90_median",
               generated=row["moisture_pathway_land_p90_median"],
               operator=">=", threshold=0.25, severity="fail",
               message="arid worlds should still have routed moisture pathways"),
    ]


def _checks_for_row(row: dict[str, Any]) -> list[dict[str, Any]]:
    preset = str(row.get("preset", "")).lower()
    label = str(row.get("label", "")).lower()
    if "earthlike" in preset or "earthlike" in label:
        return _earthlike_checks(row)
    if "waterworld" in preset or "waterworld" in label:
        return _waterworld_checks(row)
    if "arid" in preset or "arid" in label:
        return _arid_checks(row)
    return []


def _contact_sheet_arrays(summary_row: dict[str, Any]) -> dict[str, Any] | None:
    path = _array_path(summary_row)
    if path is None:
        return None
    with np.load(path, allow_pickle=False) as z:
        required = (
            "lat",
            "lon",
            "terrain__elevation_m",
            "sea_level_m",
            "atmosphere__moisture_flow_source",
            "atmosphere__moisture_flow_pathway",
            "climate__moisture_flow_network_id",
            "climate__seasonal_precipitation",
            "climate__monsoon_rainfall_corridor",
            "climate__storm_track_rainfall_corridor",
        )
        if any(key not in z for key in required):
            return None
        out = {key: np.asarray(z[key]) for key in required}
    sea = float(np.asarray(out["sea_level_m"], dtype=np.float64).ravel()[0])
    elev = np.asarray(out["terrain__elevation_m"], dtype=np.float64)
    out["land"] = elev >= sea
    return out


def _render_contact_sheet(
    summary_row: dict[str, Any],
    generated_row: dict[str, Any],
    outdir: Path,
) -> Path | None:
    arrays = _contact_sheet_arrays(summary_row)
    if arrays is None:
        return None

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from aevum.render import PRECIP_CMAP

    label = str(generated_row.get("label", "world"))
    sheets_dir = outdir / "contact_sheets"
    sheets_dir.mkdir(parents=True, exist_ok=True)
    path = sheets_dir / f"{label}_moisture_flow_network_contact_sheet.png"

    lat = np.asarray(arrays["lat"], dtype=np.float64)
    lon = np.asarray(arrays["lon"], dtype=np.float64)
    land = np.asarray(arrays["land"], dtype=bool)
    source = np.asarray(arrays["atmosphere__moisture_flow_source"], dtype=np.float64)
    pathway = np.asarray(arrays["atmosphere__moisture_flow_pathway"], dtype=np.float64)
    network_id = np.asarray(arrays["climate__moisture_flow_network_id"], dtype=np.float64)
    precip = np.asarray(arrays["climate__seasonal_precipitation"], dtype=np.float64)
    support = np.maximum(
        np.asarray(arrays["climate__monsoon_rainfall_corridor"], dtype=np.float64),
        np.asarray(arrays["climate__storm_track_rainfall_corridor"], dtype=np.float64),
    )

    seasons = ("DJF", "MAM", "JJA", "SON")
    n = max(int(lat.size), 1)
    marker_size = float(np.clip(28000.0 / n, 0.45, 10.0))
    base_color = np.where(land, "#dfdac8", "#d9edf4")

    fig, axes = plt.subplots(5, 4, figsize=(14.5, 10.5), constrained_layout=True)
    fig.suptitle(f"C4e moisture-flow-network contact sheet: {label}", fontsize=13)

    def base(ax):
        ax.scatter(lon, lat, c=base_color, s=marker_size, linewidths=0,
                   rasterized=True)
        ax.set_xlim(-180.0, 180.0)
        ax.set_ylim(-90.0, 90.0)
        ax.set_xticks([])
        ax.set_yticks([])

    def panel(ax, values, title, cmap, vmin, vmax, mask=None):
        base(ax)
        vals = np.asarray(values, dtype=np.float64)
        if mask is not None:
            vals = np.where(mask, vals, np.nan)
        finite = np.isfinite(vals)
        if finite.any():
            im = ax.scatter(lon[finite], lat[finite], c=vals[finite],
                            s=marker_size * 1.25, linewidths=0, cmap=cmap,
                            vmin=vmin, vmax=vmax, rasterized=True)
        else:
            im = ax.scatter([], [], c=[], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=8)
        return im

    precip_land = precip[:, land] if land.any() else precip
    precip_vmax = max(float(np.nanpercentile(precip_land, 98)), 100.0)
    source_vmax = max(float(np.nanpercentile(source, 98)), 0.10)
    pathway_vmax = max(float(np.nanpercentile(pathway, 98)), 0.10)
    support_vmax = max(float(np.nanpercentile(support, 98)), 0.10)
    network_positive = network_id[network_id > 0]
    network_vmax = max(float(np.nanmax(network_positive)) if network_positive.size else 1.0, 1.0)

    row_maps = [
        ("seasonal precipitation", precip, None, PRECIP_CMAP, 0.0, precip_vmax),
        ("source-ocean moisture", source, ~land, "YlGnBu", 0.0, source_vmax),
        ("land moisture pathway", pathway, land, "BuPu", 0.0, pathway_vmax),
        ("flow network id", network_id, network_id > 0.0, "tab20", 0.0, network_vmax),
        ("monsoon/storm support", support, land, "PuBuGn", 0.0, support_vmax),
    ]

    row_images = []
    for row_idx, (row_title, field, mask, cmap, vmin, vmax) in enumerate(row_maps):
        image = None
        for season_idx, season in enumerate(seasons):
            image = panel(
                axes[row_idx, season_idx],
                field[season_idx],
                f"{season} {row_title}",
                cmap,
                vmin,
                vmax,
                None if mask is None else mask[season_idx] if getattr(mask, "ndim", 1) == 2 else mask,
            )
        row_images.append(image)

    for row_idx, image in enumerate(row_images):
        fig.colorbar(image, ax=axes[row_idx, :].ravel().tolist(), shrink=0.72)

    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return path


def _render_contact_sheets(
    summary_rows: list[dict[str, Any]],
    generated_rows: list[dict[str, Any]],
    outdir: Path,
) -> list[dict[str, Any]]:
    sheets: list[dict[str, Any]] = []
    by_assets = {str(row.get("assets_dir", "")): row for row in generated_rows}
    for summary_row in summary_rows:
        generated_row = by_assets.get(str(summary_row.get("assets_dir", "")))
        if generated_row is None:
            continue
        path = _render_contact_sheet(summary_row, generated_row, outdir)
        if path is None:
            continue
        sheets.append({
            "label": str(generated_row.get("label", "")),
            "preset": str(generated_row.get("preset", "")),
            "path": str(path),
        })
    if sheets:
        manifest = {
            "schema": "aevum.earth_climate_moisture_flow_contact_sheets.v1",
            "sheet_count": int(len(sheets)),
            "sheets": sheets,
        }
        (outdir / "earth_climate_moisture_flow_contact_sheets.json").write_text(
            json.dumps(manifest, indent=2, default=_json_default),
        )
    return sheets


def _render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Earth Climate Moisture-Flow Network Gate",
        "",
        f"Schema: `{report['schema']}`",
        f"Verdict: **{report['verdict']}**",
        "",
        (
            f"Failures: `{report['failure_count']}`; warnings: "
            f"`{report['warning_count']}`; skipped: `{report['skipped_count']}`."
        ),
        "",
        "This gate evaluates C4e seasonal moisture-flow diagnostics.  It checks "
        "that Earthlike runs archive source-ocean moisture, routed land pathways, "
        "and multi-season flow-network objects; that those maps are legible "
        "corridors rather than speckles; and that waterworld/arid presets do not "
        "grow broad false-positive monsoon flow networks.",
        "",
    ]
    if report.get("contact_sheets"):
        lines.extend(["## Contact Sheets", ""])
        for sheet in report["contact_sheets"]:
            lines.append(
                f"- `{sheet['label']}` `{sheet['preset']}`: `{sheet['path']}`"
            )
        lines.append("")
    lines.extend(["## Checks", ""])
    for row in report["checks"]:
        status = "pass" if row["passed"] else row["severity"]
        lines.append(
            f"- `{status}` `{row['label']}` `{row['group']}` "
            f"`{row['metric']}` = `{row['generated']:.3f}` "
            f"{row['operator']} `{row['threshold']:.3f}`"
        )
    lines.append("")
    return "\n".join(lines)


def run_earth_climate_moisture_flow_gate(
    config: EarthClimateMoistureFlowGateConfig,
) -> dict[str, Any]:
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(Path(config.terminal_summary_json).read_text())
    summary_rows = list(summary.get("summaries", []))
    generated = [_generated_row(row) for row in summary_rows]
    contact_sheets = (
        _render_contact_sheets(summary_rows, generated, outdir)
        if config.render_contact_sheets else []
    )
    checks = [check for row in generated for check in _checks_for_row(row)]
    failures = [
        row for row in checks
        if not row["passed"] and row["severity"] == "fail"
    ]
    warnings = [
        row for row in checks
        if not row["passed"] and row["severity"] == "warn"
    ]
    skipped = [row for row in checks if row.get("skipped")]
    report = {
        "schema": SCHEMA,
        "terminal_summary_json": str(config.terminal_summary_json),
        "generated_metrics": generated,
        "checks": checks,
        "failures": failures,
        "warnings": warnings,
        "failure_count": len(failures),
        "warning_count": len(warnings),
        "skipped_count": len(skipped),
        "contact_sheets": contact_sheets,
        "contact_sheet_count": len(contact_sheets),
        "verdict": "fail" if failures else "pass",
    }
    metric_keys = sorted({
        key for row in generated for key in row.keys()
        if key not in {"kind_counts", "season_counts"}
    })
    check_keys = [
        "label", "group", "metric", "generated", "operator", "threshold",
        "severity", "passed", "skipped", "message",
    ]
    _write_csv(outdir / "earth_climate_moisture_flow_metrics.csv",
               generated, metric_keys)
    _write_csv(outdir / "earth_climate_moisture_flow_checks.csv",
               checks, check_keys)
    (outdir / "earth_climate_moisture_flow_gate_summary.json").write_text(
        json.dumps(report, indent=2, default=_json_default),
    )
    (outdir / "earth_climate_moisture_flow_gate_report.md").write_text(
        _render_markdown(report),
    )
    return report
