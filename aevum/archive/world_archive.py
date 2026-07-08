"""WorldArchive: the saved deep-time history.

We do not just keep the final map state -- we keep its *formation history*:
periodic field snapshots, the full event timeline, biological lineages, and the
ability to answer "why is this cell like this?" by assembling provenance + the
nearest causal events.
"""
from __future__ import annotations

import copy
import hashlib
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from aevum.core.events import EventBus
from aevum.core.state import WorldState


P171_REQUIRED_OBJECT_FIELDS = [
    "id",
    "kind",
    "cells",
    "birth_myr",
    "age_myr",
    "parent_process_id",
    "parent_plate_id",
    "lineage_id",
    "activity_state",
    "relief_stage",
]


@dataclass
class Frame:
    time_myr: float
    globals: dict[str, float]
    fields: dict[str, np.ndarray]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    objects: dict[str, list[dict[str, Any]]] = field(default_factory=dict)


class WorldArchive:
    DEFAULT_KEYS = [
        "terrain.elevation_m", "crust.type", "crust.age_myr", "tectonics.plate_id",
        "crust.origin", "crust.reworked_age_myr", "crust.stability",
        "crust.domain", "tectonics.continent_id", "tectonics.terrane_id",
        "tectonics.orogeny_age_myr", "terrain.province",
        "sediment.thickness_m",
        "terrain.continental_detail", "terrain.continental_detail_region_code",
        "terrain.inland_geomorphology_region_code",
        "terrain.p104f_pre_p174_lowland_support_mask",
        "terrain.p174_lowland_plain_continuity_memory",
        "terrain.p174_lowland_plain_candidate_mask",
        "terrain.p174_lowland_plain_selected_mask",
        "terrain.p174_lowland_plain_response_mask",
        "terrain.continental_province_code", "tectonics.province_parent_process",
        "terrain.old_orogen_decay_stage", "terrain.rift_margin_stage",
        "archive.wilson_cycle_phase",
        "ocean.basin_id", "ocean.margin_type", "ocean.depth_province",
        "ocean.gateway_id", "ocean.shelf_width",
        "climate.surface_temperature", "climate.precipitation", "biosphere.biome",
        "biosphere.richness", "ocean.mask",
    ]
    DEFAULT_OBJECT_KEYS = [
        "tectonics.boundary_objects",
        "tectonics.continental_provinces",
        "terrain.continental_landforms",
        "terrain.margin_landforms",
        "terrain.ocean_fabric",
        "terrain.arc_plume_landforms",
        "terrain.mountain_ranges",
        "terrain.plateau_inventory",
        "terrain.rift_margin_sequences",
    ]

    def __init__(self, world: WorldState, bus: EventBus) -> None:
        self.world = world
        self.bus = bus
        self.frames: list[Frame] = []

    # ------------------------------------------------------------------
    def capture(self, diagnostics: Optional[dict] = None,
                keys: Optional[list[str]] = None,
                object_keys: Optional[list[str]] = None) -> Frame:
        keys = keys or self.DEFAULT_KEYS
        object_keys = object_keys or self.DEFAULT_OBJECT_KEYS
        fields = {k: self.world.fields[k].copy()
                  for k in keys if k in self.world.fields}
        objects = {}
        for key in object_keys:
            if key not in self.world.objects:
                continue
            raw_objects = copy.deepcopy(self.world.objects[key])
            if not isinstance(raw_objects, list):
                continue
            objects[key] = [
                self._normalize_object_snapshot(key, obj, index)
                for index, obj in enumerate(raw_objects)
                if isinstance(obj, dict)
            ]
        frame = Frame(time_myr=self.world.time_myr,
                      globals=dict(self.world.globals),
                      fields=fields,
                      diagnostics=diagnostics or {},
                      objects=objects)
        self.frames.append(frame)
        return frame

    def _normalize_object_snapshot(
        self,
        collection: str,
        obj: dict[str, Any],
        index: int,
    ) -> dict[str, Any]:
        """Fill P171 archive-object fields without mutating live world objects."""
        out = dict(obj)
        synthesized: list[str] = []
        time_myr = float(self.world.time_myr)

        kind = str(out.get("kind") or out.get("type") or collection.rsplit(".", 1)[-1])
        if not out.get("kind"):
            out["kind"] = kind
            synthesized.append("kind")

        if not out.get("id"):
            out["id"] = self._fallback_object_id(collection, kind, out, index)
            synthesized.append("id")

        if "cells" not in out:
            cell = out.get("cell")
            out["cells"] = [int(cell)] if cell is not None else []
            synthesized.append("cells")

        birth = _first_finite(
            out.get("birth_myr"),
            out.get("formation_myr"),
            out.get("formed_myr"),
            out.get("start_myr"),
        )
        age = _first_finite(out.get("age_myr"))
        if birth is None:
            mean_age = _first_finite(out.get("mean_age_myr"))
            if age is not None:
                birth = time_myr - age
            elif mean_age is not None:
                birth = max(time_myr - mean_age, 0.0)
                age = mean_age
            else:
                birth = time_myr
        if "birth_myr" not in out:
            out["birth_myr"] = float(birth)
            synthesized.append("birth_myr")
        if age is None:
            age = max(time_myr - float(birth), 0.0)
        if "age_myr" not in out:
            out["age_myr"] = float(age)
            synthesized.append("age_myr")

        if "parent_process_id" not in out:
            parent_process = (
                out.get("parent_process_id")
                or out.get("parent_process")
                or out.get("parent_process_code")
                or _joined(out.get("parent_processes"))
                or "unknown"
            )
            out["parent_process_id"] = parent_process
            synthesized.append("parent_process_id")

        if "parent_plate_id" not in out:
            parent_plate = (
                out.get("parent_plate_id")
                if "parent_plate_id" in out
                else out.get("parent_plate_ids", out.get("plate_id"))
            )
            out["parent_plate_id"] = parent_plate
            synthesized.append("parent_plate_id")

        if "lineage_id" not in out:
            lineage = (
                out.get("lineage_id")
                if "lineage_id" in out
                else out.get("sequence_id", out.get("plateau_id", out.get("province_id", out["id"])))
            )
            out["lineage_id"] = lineage
            synthesized.append("lineage_id")

        if "activity_state" not in out:
            last_active = _first_finite(out.get("last_active_myr"))
            if last_active is not None:
                active = abs(float(last_active) - time_myr) <= 1.0e-6
                state = "active" if active else "inactive"
            elif "persistence" in out:
                state = str(out.get("persistence"))
            elif "stage" in out:
                state = str(out.get("stage"))
            else:
                state = "active"
            out["activity_state"] = state
            synthesized.append("activity_state")

        if "relief_stage" not in out:
            relief = (
                out.get("relief_stage")
                if "relief_stage" in out
                else out.get("decay_stage", out.get("mountain_class", out.get("province_class", out.get("stage", kind))))
            )
            out["relief_stage"] = relief
            synthesized.append("relief_stage")

        out["p171_required_fields_present"] = all(
            field in out for field in P171_REQUIRED_OBJECT_FIELDS
        )
        if synthesized:
            out["p171_synthesized_fields"] = sorted(set(synthesized))
        return out

    @staticmethod
    def _fallback_object_id(
        collection: str,
        kind: str,
        obj: dict[str, Any],
        index: int,
    ) -> str:
        cells = obj.get("cells", [])
        try:
            cell_values = sorted(int(c) for c in cells)[:128]
        except (TypeError, ValueError):
            cell_values = []
        if cell_values:
            payload = ",".join(str(c) for c in cell_values)
        else:
            payload = "|".join(
                str(obj.get(key, ""))
                for key in ("cell", "centroid_lat", "centroid_lon", "lat", "lon")
            ) or str(index)
        digest = hashlib.sha1(
            f"{collection}|{kind}|{payload}|{index}".encode("utf-8")
        ).hexdigest()[:12]
        return f"{collection}:{kind}:{digest}"

    # ------------------------------------------------------------------
    def timeline(self) -> list[dict]:
        return self.bus.timeline()

    def lineages(self) -> list[dict]:
        return self.world.objects.get("biosphere.lineages", [])

    def nearest_events(self, cell: int, types: Optional[list[str]] = None,
                       max_n: int = 8) -> list[dict]:
        out = []
        for e in self.bus.events:
            if types and e.type not in types:
                continue
            if e.location == cell or (isinstance(e.location, (list, np.ndarray))
                                      and cell in list(e.location)):
                out.append(e.to_dict())
        out.sort(key=lambda d: d["time_myr"])
        return out[:max_n]

    def explain_cell(self, cell: int) -> dict:
        """Assemble the causal story for one truth-layer cell."""
        w = self.world
        story: dict[str, Any] = {"cell": int(cell),
                                 "lat": round(float(w.grid.lat[cell]), 2),
                                 "lon": round(float(w.grid.lon[cell]), 2)}
        elev = float(w.get_field("terrain.elevation_m")[cell])
        story["elevation_m"] = round(elev, 1)
        story["is_ocean"] = bool(w.ocean_mask()[cell])
        story["crust_type"] = ("continental" if w.get_field("crust.type")[cell] == 1.0
                               else "oceanic")
        story["crust_age_myr"] = round(float(w.get_field("crust.age_myr")[cell]), 1)
        origin_names = {
            0: "ridge_oceanic",
            1: "primordial_continent",
            2: "arc_accreted",
            3: "collision_suture",
            4: "plume_or_impact",
            5: "craton",
        }
        origin = int(w.get_field("crust.origin", 0.0)[cell])
        story["crust_origin"] = origin_names.get(origin, f"code_{origin}")
        story["crust_reworked_age_myr"] = round(float(
            w.get_field("crust.reworked_age_myr", -1.0)[cell]), 1)
        story["crust_stability"] = round(float(
            w.get_field("crust.stability", 0.0)[cell]), 3)
        orog = float(w.get_field("tectonics.orogeny_age_myr", -1.0)[cell])
        volc = float(w.get_field("tectonics.volcanism_age_myr", -1.0)[cell])
        if orog >= 0:
            story["orogeny_formed_myr"] = round(orog, 1)
        if volc >= 0:
            story["volcanism_myr"] = round(volc, 1)
        story["tectonic_objects"] = self._cell_tectonic_objects(cell)
        story["temperature_C"] = round(float(
            w.get_field("climate.surface_temperature", 288.0)[cell]) - 273.15, 1)
        story["precip_mm_yr"] = round(float(
            w.get_field("climate.precipitation", 0.0)[cell]), 0)
        story["biome_code"] = int(w.get_field("biosphere.biome", 0.0)[cell])
        story["deposits"] = [
            {"model": d["genesis_model"], "commodities": list(d["commodity_vector"]),
             "formed_myr": d["formation_age_myr"], "grade": d["grade"],
             "burial_m": d["burial_depth_m"]}
            for d in w.objects.get("resources.deposits", []) if d["cell"] == cell]
        story["events"] = self.nearest_events(cell)
        prov = w.provenance.get("terrain.elevation_m")
        story["elevation_cause"] = prov.direct_cause if prov else ""
        return story

    def _cell_tectonic_objects(self, cell: int) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        for key in ("tectonics.boundary_objects", "tectonics.wilson_cycles",
                    "tectonics.cratons", "tectonics.lips"):
            for obj in self.world.objects.get(key, []):
                cells = obj.get("cells", [])
                if cell not in cells:
                    continue
                out.append({
                    "collection": key,
                    "id": obj.get("id"),
                    "kind": obj.get("kind"),
                    "stage": obj.get("stage"),
                })
        return out[:8]


def _first_finite(*values: Any) -> float | None:
    for value in values:
        if value is None:
            continue
        try:
            result = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(result):
            return result
    return None


def _joined(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        return "+".join(str(item) for item in value)
    return str(value) if value else ""
