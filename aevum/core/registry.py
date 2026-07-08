"""Feature Registry.

We do NOT reserve "a few hundred empty fields" on a giant World object.  Instead
every feature is registered as a *contract*: its definition, unit, spatial &
temporal representation, dependencies, producing module, valid range, fidelity
levels, uncertainty mode, conservation constraints, validation method and a
lifecycle status (reserved -> implemented -> calibrated -> validated).

This lets us "reserve all features first, fill them in one by one" while the
scheduler can already wire empty (stub) modules into a runnable pipeline.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Optional

import yaml


class Representation(str, Enum):
    """The four world-data categories (continuous field / network / object / event)
    plus convenience global scalar."""

    SCALAR_FIELD = "scalar_field"      # one value per grid cell
    VECTOR_FIELD = "vector_field"      # vector per grid cell (e.g. wind, plate velocity)
    NETWORK = "network"                # rivers, currents, plate boundaries, migration paths
    OBJECT_SET = "object_set"          # plates, volcanoes, basins, deposits, ice sheets, lineages
    EVENT_STREAM = "event_stream"      # impacts, rifts, collisions, extinctions, innovations
    GLOBAL_SCALAR = "global_scalar"    # CO2, O2, sea level, mantle temperature, luminosity


class Status(str, Enum):
    RESERVED = "reserved"
    IMPLEMENTED = "implemented"
    CALIBRATED = "calibrated"
    VALIDATED = "validated"


@dataclass
class FeatureSpec:
    feature_id: str
    representation: Representation
    domain: str                                 # e.g. surface_cells, ocean_cells, global
    unit: str
    producer: str                               # module name responsible
    dependencies: list[str] = field(default_factory=list)
    time_support: str = "snapshot"              # e.g. monthly_climatology, instantaneous, cumulative
    fidelity_levels: list[str] = field(default_factory=lambda: ["stub"])
    conservation_constraints: list[str] = field(default_factory=list)
    uncertainty: str = "none"                   # none | ensemble | interval
    valid_range: Optional[list[float]] = None   # [lo, hi] for scalar features
    validation: str = ""                        # how this feature is validated
    description: str = ""
    status: Status = Status.RESERVED

    def to_dict(self) -> dict[str, Any]:
        d = {
            "feature_id": self.feature_id,
            "representation": self.representation.value,
            "domain": self.domain,
            "unit": self.unit,
            "producer": self.producer,
            "dependencies": list(self.dependencies),
            "time_support": self.time_support,
            "fidelity_levels": list(self.fidelity_levels),
            "conservation_constraints": list(self.conservation_constraints),
            "uncertainty": self.uncertainty,
            "valid_range": self.valid_range,
            "validation": self.validation,
            "description": self.description,
            "status": self.status.value,
        }
        return d

    @staticmethod
    def from_dict(d: dict[str, Any]) -> "FeatureSpec":
        return FeatureSpec(
            feature_id=d["feature_id"],
            representation=Representation(d["representation"]),
            domain=d.get("domain", "global"),
            unit=d.get("unit", "1"),
            producer=d["producer"],
            dependencies=list(d.get("dependencies", [])),
            time_support=d.get("time_support", "snapshot"),
            fidelity_levels=list(d.get("fidelity_levels", ["stub"])),
            conservation_constraints=list(d.get("conservation_constraints", [])),
            uncertainty=d.get("uncertainty", "none"),
            valid_range=d.get("valid_range"),
            validation=d.get("validation", ""),
            description=d.get("description", ""),
            status=Status(d.get("status", "reserved")),
        )


class FeatureRegistry:
    def __init__(self) -> None:
        self._specs: dict[str, FeatureSpec] = {}

    # -- registration ----------------------------------------------------
    def register(self, spec: FeatureSpec) -> FeatureSpec:
        if spec.feature_id in self._specs:
            raise ValueError(f"feature already registered: {spec.feature_id}")
        self._specs[spec.feature_id] = spec
        return spec

    def register_many(self, specs: Iterable[FeatureSpec]) -> None:
        for s in specs:
            self.register(s)

    # -- access ----------------------------------------------------------
    def get(self, feature_id: str) -> FeatureSpec:
        return self._specs[feature_id]

    def __contains__(self, feature_id: str) -> bool:
        return feature_id in self._specs

    def __iter__(self):
        return iter(self._specs.values())

    def __len__(self) -> int:
        return len(self._specs)

    def ids(self) -> list[str]:
        return list(self._specs)

    def by_producer(self, producer: str) -> list[FeatureSpec]:
        return [s for s in self._specs.values() if s.producer == producer]

    def by_status(self, status: Status) -> list[FeatureSpec]:
        return [s for s in self._specs.values() if s.status == status]

    # -- io --------------------------------------------------------------
    def load_yaml_dir(self, directory: str | Path) -> int:
        directory = Path(directory)
        count = 0
        for path in sorted(directory.glob("*.y*ml")):
            with open(path, "r", encoding="utf-8") as fh:
                docs = yaml.safe_load(fh)
            entries = docs if isinstance(docs, list) else docs.get("features", [])
            for entry in entries:
                self.register(FeatureSpec.from_dict(entry))
                count += 1
        return count

    def dump_yaml(self, path: str | Path) -> None:
        data = {"features": [s.to_dict() for s in self._specs.values()]}
        with open(path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(data, fh, sort_keys=False, allow_unicode=True)

    # -- validation ------------------------------------------------------
    def validate(self) -> list[str]:
        """Return *fatal* contract problems (empty list == healthy registry).

        Dependency cycles are NOT fatal: they are physical feedback loops
        (climate<->weathering<->CO2, terrain<->precipitation<->erosion,
        nutrients<->biomass<->burial) that the scheduler resolves by iteration
        within each macro time slice.  Use :meth:`feedback_loops` to list them.
        """
        problems: list[str] = []
        known = set(self._specs)
        for spec in self._specs.values():
            for dep in spec.dependencies:
                if dep not in known:
                    problems.append(f"{spec.feature_id}: unknown dependency '{dep}'")
        return problems

    def feedback_loops(self) -> list[str]:
        """Return detected dependency cycles (expected, informational)."""
        return self._detect_cycles()

    def _detect_cycles(self) -> list[str]:
        WHITE, GREY, BLACK = 0, 1, 2
        color = {fid: WHITE for fid in self._specs}
        problems: list[str] = []

        def visit(fid: str, stack: list[str]) -> None:
            color[fid] = GREY
            spec = self._specs.get(fid)
            if spec is not None:
                for dep in spec.dependencies:
                    if dep not in self._specs:
                        continue
                    if color[dep] == GREY:
                        cyc = " -> ".join(stack + [fid, dep])
                        problems.append(f"dependency cycle: {cyc}")
                    elif color[dep] == WHITE:
                        visit(dep, stack + [fid])
            color[fid] = BLACK

        for fid in self._specs:
            if color[fid] == WHITE:
                visit(fid, [])
        return problems

    def status_summary(self) -> dict[str, int]:
        out = {s.value: 0 for s in Status}
        for spec in self._specs.values():
            out[spec.status.value] += 1
        return out
