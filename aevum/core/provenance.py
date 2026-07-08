"""Provenance: every output can answer "why is it like this?".

For each feature we record which module produced it, at what fidelity, with what
uncertainty, its direct cause (a short human string) and the upstream event ids.
The map compiler exposes this so a player clicking a tile gets a causal story.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class Provenance:
    feature_id: str
    producer: str
    fidelity: str = "stub"
    unit: str = "1"
    uncertainty: Optional[tuple[float, float]] = None
    direct_cause: str = ""
    upstream_events: list[int] = field(default_factory=list)
    updated_at_myr: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "feature_id": self.feature_id,
            "producer": self.producer,
            "fidelity": self.fidelity,
            "unit": self.unit,
            "uncertainty": list(self.uncertainty) if self.uncertainty else None,
            "direct_cause": self.direct_cause,
            "upstream_events": list(self.upstream_events),
            "updated_at_myr": round(self.updated_at_myr, 4),
        }


class ProvenanceStore:
    def __init__(self) -> None:
        self._by_feature: dict[str, Provenance] = {}

    def record(self, prov: Provenance) -> None:
        self._by_feature[prov.feature_id] = prov

    def get(self, feature_id: str) -> Optional[Provenance]:
        return self._by_feature.get(feature_id)

    def all(self) -> dict[str, Provenance]:
        return dict(self._by_feature)
