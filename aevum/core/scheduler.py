"""Deep-time scheduler: multi-rate, adaptive, event-aware.

We never integrate 4.5 Gyr year-by-year.  Modules run at their own cadences;
expensive modules (climate) re-solve only when their drivers (land/sea
distribution, orbit, atmosphere, topography) drift past a threshold -- otherwise
the previous solution is reused.  The macro time-step itself is adaptive: fine
early (hot, tectonically vigorous) and coarse later.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import time
from typing import Any, Callable, Optional

import numpy as np

from aevum.core.events import Event, EventBus
from aevum.core.module import Module, StepResult
from aevum.core.rng import RNGKey
from aevum.core.state import WorldState


@dataclass
class ScheduledModule:
    module: Module
    interval_myr: float
    # Optional re-solve trigger: returns True if the module must run now even if
    # not yet due (e.g. climate when geography changed a lot).
    trigger: Optional[Callable[[WorldState], bool]] = None
    last_run_myr: float = -np.inf
    runs: int = 0


@dataclass
class StepRecord:
    time_myr: float
    dt_myr: float
    ran: list[str]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    module_seconds: dict[str, float] = field(default_factory=dict)


class DeepTimeScheduler:
    def __init__(self, world: WorldState, bus: EventBus,
                 t_end_myr: float,
                 dt_start_myr: float = 10.0,
                 dt_end_myr: float = 50.0) -> None:
        self.world = world
        self.bus = bus
        self.t_end = t_end_myr
        self.dt_start = dt_start_myr
        self.dt_end = dt_end_myr
        self.modules: list[ScheduledModule] = []
        self.history: list[StepRecord] = []
        self._event_counter = 0

    def add(self, module: Module, interval_myr: Optional[float] = None,
            trigger: Optional[Callable[[WorldState], bool]] = None) -> None:
        self.modules.append(ScheduledModule(
            module=module,
            interval_myr=interval_myr if interval_myr is not None else module.interval_myr,
            trigger=trigger,
        ))

    # -- adaptive macro step --------------------------------------------
    def _dt(self, t: float) -> float:
        frac = min(1.0, t / max(self.t_end, 1e-9))
        return float(self.dt_start + (self.dt_end - self.dt_start) * frac)

    # -- lifecycle -------------------------------------------------------
    def init(self) -> None:
        for sm in self.modules:
            key = RNGKey(self.world.spec.seed, sm.module.name, 0.0)
            sm.module.init_state(self.world, key)

    def run(self, on_step: Optional[Callable[[StepRecord], None]] = None) -> None:
        self.init()
        t = 0.0
        while t < self.t_end - 1e-9:
            dt = min(self._dt(t), self.t_end - t)
            t_next = t + dt
            self.world.time_myr = t_next
            ran: list[str] = []
            diag: dict[str, Any] = {}
            module_seconds: dict[str, float] = {}
            final_step = t_next >= self.t_end - 1e-9
            for sm in self.modules:
                due = (
                    final_step
                    or (t_next - sm.last_run_myr) >= sm.interval_myr - 1e-9
                )
                triggered = sm.trigger(self.world) if sm.trigger else False
                if not (due or triggered):
                    continue
                eff_dt = t_next - max(sm.last_run_myr, 0.0)
                self._event_counter += 1
                key = RNGKey(self.world.spec.seed, sm.module.name, t_next,
                             self._event_counter)
                forcing = {"due": due, "triggered": triggered}
                module_start = time.perf_counter()
                result = sm.module.step(self.world, t_next, eff_dt, forcing, key)
                module_seconds[sm.module.name] = float(
                    time.perf_counter() - module_start
                )
                self._apply(result, sm.module, t_next)
                sm.last_run_myr = t_next
                sm.runs += 1
                ran.append(sm.module.name)
                if result.diagnostics:
                    diag[sm.module.name] = result.diagnostics
            rec = StepRecord(
                time_myr=t_next,
                dt_myr=dt,
                ran=ran,
                diagnostics=diag,
                module_seconds=module_seconds,
            )
            self.history.append(rec)
            if on_step is not None:
                on_step(rec)
            t = t_next

    # -- delta application ----------------------------------------------
    def _apply(self, result: StepResult, module: Module, t: float) -> None:
        for event in result.events:
            self.bus.emit(event)
        delta = result.state_delta
        for key, value in delta.get("fields", {}).items():
            self.world.fields[key] = value
        for key, value in delta.get("globals", {}).items():
            self.world.globals[key] = value
        for key, value in delta.get("networks", {}).items():
            self.world.networks[key] = value
        for key, value in delta.get("objects", {}).items():
            self.world.objects[key] = value
