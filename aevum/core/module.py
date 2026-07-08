"""Unified module interface.

Every process module implements the same contract:

    step(world, t, dt, forcing, rng_key) -> StepResult(state_delta, events, diagnostics)

The engine applies ``state_delta`` to the world, appends ``events`` to the bus and
records ``diagnostics``.  Modules declare which features they ``produces`` and
``consumes`` (must match the registry) plus a base ``interval_myr`` cadence and a
``fidelity`` level.  A stub module simply returns an empty result, which keeps the
whole pipeline runnable while features are still ``reserved``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from aevum.core.events import Event
from aevum.core.rng import RNGKey
from aevum.core.state import WorldState


@dataclass
class StepResult:
    state_delta: dict[str, Any] = field(default_factory=dict)
    events: list[Event] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)


class Module:
    name: str = "module"
    produces: list[str] = []
    consumes: list[str] = []
    fidelity: str = "stub"
    interval_myr: float = 100.0

    def init_state(self, world: WorldState, rng_key: RNGKey) -> None:
        """Optional: seed initial fields/objects at t=0."""

    def step(self, world: WorldState, t: float, dt: float,
             forcing: dict[str, Any], rng_key: RNGKey) -> StepResult:
        raise NotImplementedError

    # -- convenience -----------------------------------------------------
    def key(self, world: WorldState, t: float, event_index: int = 0) -> RNGKey:
        return RNGKey(world.spec.seed, self.name, t, event_index)


class StubModule(Module):
    """A do-nothing module used to reserve a slot in the pipeline."""

    def __init__(self, name: str, produces: list[str] | None = None,
                 interval_myr: float = 100.0):
        self.name = name
        self.produces = produces or []
        self.interval_myr = interval_myr

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        return StepResult()
