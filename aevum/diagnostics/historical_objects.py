"""P171 historical geomorphology object persistence diagnostics."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from aevum.archive.world_archive import (
    P171_REQUIRED_OBJECT_FIELDS,
    WorldArchive,
)


SCHEMA = "aevum.p171_historical_object_persistence.v1"


def historical_object_persistence_summary(
    world: Any,
    archive: Any,
    *,
    object_keys: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Summarize archived process object completeness and ID persistence."""
    keys = tuple(object_keys or tuple(WorldArchive.DEFAULT_OBJECT_KEYS))
    frames = list(getattr(archive, "frames", []) or [])
    rows: list[dict[str, Any]] = []
    missing_field_counts: dict[str, int] = {}
    collection_counts: dict[str, int] = {}
    observations_by_key: dict[tuple[str, str], list[float]] = {}
    total_objects = 0
    total_required_slots = 0
    missing_required_slots = 0

    for frame in frames:
        frame_time = float(getattr(frame, "time_myr", 0.0))
        objects_by_key = getattr(frame, "objects", {}) or {}
        objects_by_key = objects_by_key if isinstance(objects_by_key, dict) else {}
        for collection in keys:
            objects = objects_by_key.get(collection, [])
            if not isinstance(objects, list):
                objects = []
            object_count = 0
            missing_in_collection: dict[str, int] = {}
            kind_counts: dict[str, int] = {}
            for obj in objects:
                if not isinstance(obj, dict):
                    continue
                object_count += 1
                total_objects += 1
                kind = str(obj.get("kind", "unknown") or "unknown")
                kind_counts[kind] = kind_counts.get(kind, 0) + 1
                obj_id = str(obj.get("id", ""))
                if obj_id:
                    observations_by_key.setdefault((collection, obj_id), []).append(frame_time)
                for field in P171_REQUIRED_OBJECT_FIELDS:
                    total_required_slots += 1
                    if field not in obj:
                        missing_required_slots += 1
                        key = f"{collection}.{field}"
                        missing_field_counts[key] = missing_field_counts.get(key, 0) + 1
                        missing_in_collection[field] = missing_in_collection.get(field, 0) + 1
            collection_counts[collection] = collection_counts.get(collection, 0) + object_count
            if object_count or collection in objects_by_key:
                rows.append({
                    "time_myr": frame_time,
                    "collection": collection,
                    "object_count": int(object_count),
                    "kind_counts": dict(sorted(kind_counts.items())),
                    "missing_required_field_counts": dict(sorted(missing_in_collection.items())),
                    "required_fields_complete": bool(not missing_in_collection),
                })

    recurring: list[dict[str, Any]] = []
    for (collection, obj_id), times in observations_by_key.items():
        unique_times = sorted({float(t) for t in times})
        if len(unique_times) <= 1:
            continue
        recurring.append({
            "collection": collection,
            "id": obj_id,
            "frame_count": int(len(unique_times)),
            "first_time_myr": float(unique_times[0]),
            "last_time_myr": float(unique_times[-1]),
        })
    recurring.sort(key=lambda item: (
        str(item["collection"]),
        str(item["id"]),
        float(item["first_time_myr"]),
    ))

    unique_object_ids = len(observations_by_key)
    required_fields_complete = missing_required_slots == 0
    return {
        "schema": SCHEMA,
        "context": {
            "spec_name": str(getattr(getattr(world, "spec", None), "name", "")),
            "seed": int(getattr(getattr(world, "spec", None), "seed", 0)),
            "final_time_myr": float(getattr(world, "time_myr", 0.0)),
        },
        "required_object_fields": list(P171_REQUIRED_OBJECT_FIELDS),
        "object_keys": list(keys),
        "frame_count": int(len(frames)),
        "object_collection_row_count": int(len(rows)),
        "total_object_observations": int(total_objects),
        "unique_object_id_count": int(unique_object_ids),
        "recurring_object_id_count": int(len(recurring)),
        "required_field_slot_count": int(total_required_slots),
        "missing_required_field_slot_count": int(missing_required_slots),
        "required_fields_complete": bool(required_fields_complete),
        "missing_required_field_counts": dict(sorted(missing_field_counts.items())),
        "collection_object_counts": dict(sorted(collection_counts.items())),
        "recurring_objects_sample": recurring[:32],
        "frame_collection_rows": rows,
        "acceptance": {
            "archive_frames_present": bool(len(frames) > 0),
            "archive_objects_present": bool(total_objects > 0),
            "required_fields_complete": bool(required_fields_complete),
            "persistence_checked": bool(len(frames) >= 2),
            "recurring_object_ids_present": bool(len(recurring) > 0),
            "generation_behavior_changed": False,
        },
    }


def write_historical_object_audit(
    world: Any,
    archive: Any,
    outdir: str | Path,
    *,
    filename: str = "p171_historical_object_persistence_audit.json",
) -> dict[str, Any]:
    """Write the P171 object persistence audit JSON and return it."""
    summary = historical_object_persistence_summary(world, archive)
    path = Path(outdir)
    path.mkdir(parents=True, exist_ok=True)
    (path / filename).write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary
