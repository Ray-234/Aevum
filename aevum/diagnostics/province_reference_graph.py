"""Derived real-Earth province reference graph for tectonics benchmarks.

The graph is a compact fixture extracted from the P73 case-study sketches.  It
does not bundle raw GIS vectors; it makes the province class, process, and
adjacency expectations executable while raw-vector extraction remains a later
stage.
"""
from __future__ import annotations

from collections import Counter, defaultdict, deque
import hashlib
import json
from typing import Any

from aevum.diagnostics.earth_case_studies import (
    REQUIRED_FEATURE_CLASSES,
    REQUIRED_PROCESSES,
    earth_case_study_calibration_summary,
)


SCHEMA = "aevum.province_reference_graph.v1"

REQUIRED_CLASS_EDGES = {
    ("active_orogen", "foreland_basin"),
    ("active_orogen", "old_orogen"),
    ("intracratonic_basin", "platform"),
    ("intracratonic_basin", "shield"),
    ("old_orogen", "passive_margin_lowland"),
    ("old_orogen", "platform"),
    ("passive_margin_lowland", "shield"),
    ("platform", "shield"),
    ("rift_system", "volcanic_lip_plateau"),
}


def _pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((str(a), str(b))))


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _component_count(nodes: set[str], edges: tuple[tuple[str, str], ...]) -> int:
    if not nodes:
        return 0
    adjacency: dict[str, set[str]] = {node: set() for node in nodes}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    seen: set[str] = set()
    count = 0
    for node in sorted(nodes):
        if node in seen:
            continue
        count += 1
        queue: deque[str] = deque([node])
        seen.add(node)
        while queue:
            current = queue.popleft()
            for nxt in adjacency.get(current, ()):
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return count


def _case_reference_graph(case: dict[str, Any]) -> dict[str, Any]:
    province_by_id = {
        str(province["id"]): province for province in case["provinces"]
    }
    node_ids = set(province_by_id)
    raw_edges = tuple(
        (str(edge[0]), str(edge[1])) for edge in case["adjacency"]
    )
    undirected_edges = tuple(_pair(a, b) for a, b in raw_edges)
    edge_counts = Counter(undirected_edges)
    duplicate_edges = tuple(sorted(edge for edge, count in edge_counts.items() if count > 1))
    invalid_edges = tuple(
        edge for edge in undirected_edges
        if edge[0] not in node_ids or edge[1] not in node_ids
    )
    incident = {node for edge in undirected_edges for node in edge}
    isolated_nodes = tuple(sorted(node_ids - incident))
    classes = Counter(str(province["class"]) for province in province_by_id.values())
    processes = sorted({
        str(process)
        for province in province_by_id.values()
        for process in province["parent_processes"]
    })
    class_edges = Counter(
        _pair(province_by_id[a]["class"], province_by_id[b]["class"])
        for a, b in undirected_edges
        if a in province_by_id and b in province_by_id
    )
    node_count = len(node_ids)
    edge_count = len(undirected_edges)
    return {
        "case_id": str(case["id"]),
        "label": str(case["label"]),
        "source_ids": tuple(case["source_ids"]),
        "node_count": node_count,
        "edge_count": edge_count,
        "component_count": _component_count(node_ids, undirected_edges),
        "isolated_node_count": len(isolated_nodes),
        "isolated_nodes": isolated_nodes,
        "duplicate_edge_count": len(duplicate_edges),
        "duplicate_edges": duplicate_edges,
        "invalid_edge_count": len(invalid_edges),
        "invalid_edges": invalid_edges,
        "class_count": len(classes),
        "classes": tuple(sorted(classes)),
        "class_counts": dict(sorted(classes.items())),
        "largest_class_fraction": (
            float(max(classes.values()) / node_count) if node_count else 0.0
        ),
        "parent_process_count": len(processes),
        "parent_processes": tuple(processes),
        "class_edge_count": len(class_edges),
        "class_edges": {
            "|".join(edge): int(count) for edge, count in sorted(class_edges.items())
        },
        "edge_density": (
            float(2.0 * edge_count / (node_count * (node_count - 1)))
            if node_count > 1 else 0.0
        ),
        "acceptance": {
            "connected": _component_count(node_ids, undirected_edges) == 1,
            "no_isolated_nodes": not isolated_nodes,
            "no_duplicate_edges": not duplicate_edges,
            "edges_valid": not invalid_edges,
            "multi_class_graph": len(classes) >= 4,
            "source_ids_known": not case["missing_source_ids"],
            "metric_envelope_valid": (
                not case["missing_envelope_keys"]
                and not case["invalid_envelope_keys"]
            ),
        },
    }


def province_reference_graph_summary() -> dict[str, Any]:
    """Return a compact, deterministic real-Earth province graph fixture."""
    case_summary = earth_case_study_calibration_summary()
    cases = tuple(_case_reference_graph(case) for case in case_summary["cases"])

    all_classes = sorted({
        province_class
        for case in cases
        for province_class in case["classes"]
    })
    all_processes = sorted({
        process
        for case in cases
        for process in case["parent_processes"]
    })
    class_edge_counts: Counter[tuple[str, str]] = Counter()
    source_ids = set()
    for case in cases:
        source_ids.update(str(source_id) for source_id in case["source_ids"])
        for edge, count in case["class_edges"].items():
            a, b = edge.split("|", maxsplit=1)
            class_edge_counts[(a, b)] += int(count)

    missing_feature_classes = tuple(sorted(set(REQUIRED_FEATURE_CLASSES) - set(all_classes)))
    missing_parent_processes = tuple(sorted(set(REQUIRED_PROCESSES) - set(all_processes)))
    missing_class_edges = tuple(sorted(REQUIRED_CLASS_EDGES - set(class_edge_counts)))
    acceptance = {
        "case_study_calibration_ready": (
            case_summary["status"] == "case_study_calibration_ready"
        ),
        "small_derived_fixture": True,
        "raw_vectors_not_stored": True,
        "direct_vector_extraction_marked_pending": True,
        "feature_class_not_exact_geography_policy": True,
        "required_feature_classes_covered": not missing_feature_classes,
        "required_parent_processes_covered": not missing_parent_processes,
        "required_class_edges_covered": not missing_class_edges,
        "all_case_sources_known": all(
            case["acceptance"]["source_ids_known"] for case in case_summary["cases"]
        ),
        "all_case_graphs_connected": all(
            case["acceptance"]["connected"] for case in cases
        ),
        "all_case_edges_valid": all(
            case["acceptance"]["edges_valid"] for case in cases
        ),
        "no_duplicate_case_edges": all(
            case["acceptance"]["no_duplicate_edges"] for case in cases
        ),
        "no_isolated_provinces": all(
            case["acceptance"]["no_isolated_nodes"] for case in cases
        ),
    }
    summary = {
        "schema": SCHEMA,
        "source_case_study_schema": case_summary["schema"],
        "source_case_study_status": case_summary["status"],
        "case_count": len(cases),
        "node_count": sum(case["node_count"] for case in cases),
        "edge_count": sum(case["edge_count"] for case in cases),
        "class_count": len(all_classes),
        "parent_process_count": len(all_processes),
        "source_id_count": len(source_ids),
        "class_edge_count": len(class_edge_counts),
        "classes": tuple(all_classes),
        "parent_processes": tuple(all_processes),
        "source_ids": tuple(sorted(source_ids)),
        "class_edge_counts": {
            "|".join(edge): int(count)
            for edge, count in sorted(class_edge_counts.items())
        },
        "required_class_edges": tuple(sorted(REQUIRED_CLASS_EDGES)),
        "missing_required_class_edges": missing_class_edges,
        "missing_required_feature_classes": missing_feature_classes,
        "missing_required_parent_processes": missing_parent_processes,
        "cases": cases,
        "extraction_policy": {
            "derived_from": "aevum.diagnostics.earth_case_studies.CASE_STUDIES",
            "raw_vector_extraction": False,
            "raw_vector_extraction_pending": True,
            "exact_geography_required": False,
            "feature_class_required": True,
            "intended_use": (
                "Regression target for generated province graphs and future "
                "raw-GIS extraction."
            ),
        },
        "acceptance": acceptance,
        "next_gates": (
            "P80.generated_province_graph_reference_comparison",
            "P81.real_earth_landform_catalog_expansion",
        ),
    }
    summary["status"] = (
        "province_reference_graph_ready"
        if all(acceptance.values())
        else "province_reference_graph_incomplete"
    )
    summary["fixture_digest"] = _digest(summary)
    return summary
