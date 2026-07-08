"""Offline source-ledger schema for real-Earth reference extraction.

The P76 ledger deliberately does not download reference data.  It makes the
metadata contract explicit before any large rasters or vectors are acquired.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

from aevum.diagnostics.physiographic_reference import PHASES, REFERENCE_SOURCES


SCHEMA = "aevum.reference_source_ledger.v1"
LEDGER_VERSION = "2026-06-27.p76"

REQUIRED_LEDGER_FIELDS: tuple[str, ...] = (
    "source_id",
    "category",
    "phase_ids",
    "source_url",
    "source_kind",
    "source_version_note",
    "license_status",
    "acquisition_status",
    "extraction_status",
    "raw_data_policy",
    "local_storage_policy",
    "projection_resolution_note",
    "checksum_status",
    "derived_metric_targets",
)

RAW_DATA_POLICY = (
    "Do not commit raw external rasters or vectors.  Keep raw acquisitions "
    "outside git and store only reproducible scripts plus small derived metric "
    "fixtures in the repo."
)


def _ledger_entry(source: dict[str, Any]) -> dict[str, Any]:
    source_id = str(source["id"])
    url = str(source["url"])
    uses = tuple(str(use) for use in source.get("uses", ()))
    phases = tuple(str(phase) for phase in source.get("phases", ()))
    is_internal = url.startswith("repo://")
    if is_internal:
        return {
            "source_id": source_id,
            "category": str(source["category"]),
            "phase_ids": phases,
            "source_url": url,
            "source_kind": "internal_generated_archive",
            "source_version_note": "archived generated baseline: 20260626_c",
            "license_status": "internal_project_output",
            "acquisition_status": "available_local_archive",
            "extraction_status": "summary_metrics_available",
            "raw_data_policy": "Internal generated artifacts may be archived; no external raw reference data is bundled.",
            "local_storage_policy": "Use archived summary JSON and rendered PNG assets under the referenced repo path.",
            "projection_resolution_note": "Aevum spherical grid, 8000-cell generated audit; rendered diagnostic PNGs are equirectangular.",
            "checksum_status": "archive_path_and_summary_digest_recorded_by_benchmark",
            "derived_metric_targets": uses,
        }
    return {
        "source_id": source_id,
        "category": str(source["category"]),
        "phase_ids": phases,
        "source_url": url,
        "source_kind": "external_reference",
        "source_version_note": "record exact source release/version during Stage A extraction",
        "license_status": "license_review_required_before_download",
        "acquisition_status": "not_downloaded_by_default",
        "extraction_status": "pending_stage_a_extraction",
        "raw_data_policy": RAW_DATA_POLICY,
        "local_storage_policy": "Store raw data outside git; commit only source ledger, extraction code, and small derived metric JSON fixtures.",
        "projection_resolution_note": "record native projection, datum, vertical reference, horizontal resolution, and any downsampling during extraction",
        "checksum_status": "pending_until_acquisition",
        "derived_metric_targets": uses,
    }


def source_ledger_summary() -> dict[str, Any]:
    entries = tuple(_ledger_entry(source) for source in REFERENCE_SOURCES)
    source_ids = [str(entry["source_id"]) for entry in entries]
    p70_source_ids = [str(source["id"]) for source in REFERENCE_SOURCES]
    phase_ids = {str(phase["id"]) for phase in PHASES}
    category_counts = Counter(str(entry["category"]) for entry in entries)
    missing_fields = {
        str(entry["source_id"]): tuple(
            field
            for field in REQUIRED_LEDGER_FIELDS
            if field not in entry or entry[field] in ("", (), None)
        )
        for entry in entries
    }
    missing_fields = {
        source_id: fields
        for source_id, fields in missing_fields.items()
        if fields
    }
    invalid_phase_refs = {
        str(entry["source_id"]): tuple(
            phase for phase in entry["phase_ids"] if phase not in phase_ids
        )
        for entry in entries
    }
    invalid_phase_refs = {
        source_id: phases
        for source_id, phases in invalid_phase_refs.items()
        if phases
    }
    covered_phases = {phase for entry in entries for phase in entry["phase_ids"]}
    acceptance = {
        "all_required_fields_present": not missing_fields,
        "source_ids_match_p70_inventory": set(source_ids) == set(p70_source_ids),
        "unique_source_ids": len(set(source_ids)) == len(source_ids),
        "all_declared_phases_valid": not invalid_phase_refs,
        "all_p70_phases_covered": phase_ids.issubset(covered_phases),
        "license_status_explicit": all(bool(entry["license_status"]) for entry in entries),
        "source_version_note_explicit": all(bool(entry["source_version_note"]) for entry in entries),
        "projection_resolution_note_explicit": all(
            bool(entry["projection_resolution_note"]) for entry in entries),
        "raw_data_policy_explicit": all(bool(entry["raw_data_policy"]) for entry in entries),
        "acquisition_separated_from_benchmarks": all(
            str(entry["acquisition_status"]) in {
                "available_local_archive",
                "not_downloaded_by_default",
            }
            for entry in entries
        ),
        "derived_metric_targets_present": all(
            bool(entry["derived_metric_targets"]) for entry in entries),
        "no_external_raw_data_required": all(
            entry["source_kind"] == "internal_generated_archive"
            or entry["acquisition_status"] == "not_downloaded_by_default"
            for entry in entries
        ),
    }
    ready = all(acceptance.values())
    return {
        "schema": SCHEMA,
        "ledger_version": LEDGER_VERSION,
        "status": "source_ledger_schema_ready" if ready else "source_ledger_schema_incomplete",
        "required_fields": REQUIRED_LEDGER_FIELDS,
        "raw_data_policy": RAW_DATA_POLICY,
        "source_count": len(entries),
        "external_source_count": int(sum(
            entry["source_kind"] == "external_reference" for entry in entries)),
        "internal_source_count": int(sum(
            entry["source_kind"] == "internal_generated_archive" for entry in entries)),
        "category_count": len(category_counts),
        "category_counts": dict(sorted(category_counts.items())),
        "missing_required_fields": missing_fields,
        "invalid_phase_refs": invalid_phase_refs,
        "duplicate_source_ids": sorted(
            source_id for source_id, count in Counter(source_ids).items() if count > 1),
        "entries": entries,
        "acceptance": acceptance,
    }
