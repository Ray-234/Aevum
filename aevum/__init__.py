"""aevum: planetary deep-time evolution engine.

Three strictly separated layers:

1. Truth layer   -- physically consistent planet state (``WorldState``).
2. History layer -- events, lineages, strata, deposits (``WorldArchive``).
3. Game layer    -- compiled strategy hex map (``MapCompiler``).

Game balancing may only happen in layer 3 and must never silently mutate the
truth layer.  The internal timeline is *time since planet formation* (Myr);
geological era names are only display labels of the Earth preset scenario.
"""

__version__ = "0.1.0"

from aevum.core.units import Dimension, Conserved, CONSTANTS
from aevum.core.registry import FeatureRegistry, FeatureSpec, Representation, Status
from aevum.core.state import WorldState
from aevum.core.grid import SphereGrid
from aevum.core.module import Module, StepResult
from aevum.core.scheduler import DeepTimeScheduler, ScheduledModule
from aevum.core.events import Event, EventBus
from aevum.spec.planet_spec import PlanetSpec

__all__ = [
    "__version__",
    "Dimension",
    "Conserved",
    "CONSTANTS",
    "FeatureRegistry",
    "FeatureSpec",
    "Representation",
    "Status",
    "WorldState",
    "SphereGrid",
    "Module",
    "StepResult",
    "DeepTimeScheduler",
    "ScheduledModule",
    "Event",
    "EventBus",
    "PlanetSpec",
]
