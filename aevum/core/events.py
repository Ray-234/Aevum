"""Event format and event bus (the catastrophe / history backbone).

Events are first-class citizens of the history layer.  Each event has a stable
id, a type, a time, an optional spatial location, a magnitude, free-form params,
the module that emitted it, and a list of *cause* event ids so the archive can
reconstruct causal chains (e.g. ``rift -> ocean_basin -> collision -> orogeny``).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Any, Optional

_event_counter = itertools.count(1)


@dataclass
class Event:
    type: str
    time_myr: float
    producer: str
    location: Optional[Any] = None      # cell id, list of cells, region label, or None=global
    magnitude: float = 0.0
    params: dict[str, Any] = field(default_factory=dict)
    causes: list[int] = field(default_factory=list)
    id: int = field(default_factory=lambda: next(_event_counter))

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "type": self.type,
            "time_myr": round(self.time_myr, 4),
            "producer": self.producer,
            "location": self.location,
            "magnitude": self.magnitude,
            "params": self.params,
            "causes": self.causes,
        }


class EventBus:
    """Collects emitted events for the current run."""

    def __init__(self) -> None:
        self._events: list[Event] = []

    def emit(self, event: Event) -> Event:
        self._events.append(event)
        return event

    def extend(self, events: list[Event]) -> None:
        self._events.extend(events)

    @property
    def events(self) -> list[Event]:
        return self._events

    def by_type(self, type_: str) -> list[Event]:
        return [e for e in self._events if e.type == type_]

    def in_window(self, t0: float, t1: float) -> list[Event]:
        return [e for e in self._events if t0 <= e.time_myr <= t1]

    def timeline(self) -> list[dict[str, Any]]:
        return [e.to_dict() for e in sorted(self._events, key=lambda e: e.time_myr)]
