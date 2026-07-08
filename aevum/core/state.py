"""WorldState: the truth-layer container.

World data is split into the four categories from the project plan instead of one
giant flat object:

  * ``fields``   -- continuous scalar/vector fields over grid cells
  * ``networks`` -- rivers, currents, plate boundaries, migration paths
  * ``objects``  -- plates, volcanoes, basins, deposits, ice sheets, lineages
  * ``globals``  -- single global scalars (CO2, O2, sea level, mantle temperature)

Events live in the :class:`~aevum.core.events.EventBus` (history layer) but are
referenced from provenance.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.provenance import Provenance, ProvenanceStore


@dataclass
class WorldState:
    grid: SphereGrid
    spec: Any                                   # PlanetSpec (avoid import cycle)
    time_myr: float = 0.0

    fields: dict[str, np.ndarray] = field(default_factory=dict)
    networks: dict[str, Any] = field(default_factory=dict)
    objects: dict[str, list] = field(default_factory=dict)
    globals: dict[str, float] = field(default_factory=dict)

    provenance: ProvenanceStore = field(default_factory=ProvenanceStore)

    # -- field helpers ---------------------------------------------------
    @property
    def n_cells(self) -> int:
        return self.grid.n

    def field(self, name: str) -> np.ndarray:
        return self.fields[name]

    def get_field(self, name: str, default: float = 0.0) -> np.ndarray:
        if name not in self.fields:
            return np.full(self.n_cells, default, dtype=np.float64)
        return self.fields[name]

    def set_field(self, name: str, values: np.ndarray,
                  provenance: Optional[Provenance] = None) -> None:
        arr = np.asarray(values, dtype=np.float64)
        self.fields[name] = arr
        if provenance is not None:
            provenance.updated_at_myr = self.time_myr
            self.provenance.record(provenance)

    def ensure_field(self, name: str, fill: float = 0.0) -> np.ndarray:
        if name not in self.fields:
            self.fields[name] = np.full(self.n_cells, fill, dtype=np.float64)
        return self.fields[name]

    # -- object helpers --------------------------------------------------
    def object_set(self, name: str) -> list:
        return self.objects.setdefault(name, [])

    # -- global helpers --------------------------------------------------
    def g(self, name: str, default: float = 0.0) -> float:
        return self.globals.get(name, default)

    def set_g(self, name: str, value: float) -> None:
        self.globals[name] = float(value)

    # -- derived masks ---------------------------------------------------
    @property
    def sea_level(self) -> float:
        return self.g("ocean.sea_level_m", 0.0)

    def ocean_mask(self) -> np.ndarray:
        """Boolean: cell surface is below sea level."""
        elev = self.get_field("terrain.elevation_m")
        return elev < self.sea_level

    def land_mask(self) -> np.ndarray:
        return ~self.ocean_mask()

    def water_depth(self) -> np.ndarray:
        elev = self.get_field("terrain.elevation_m")
        return np.maximum(self.sea_level - elev, 0.0)

    def land_fraction(self) -> float:
        land = self.land_mask()
        area = self.grid.cell_area
        return float(area[land].sum() / area.sum())

    def snapshot_fields(self) -> dict[str, np.ndarray]:
        return {k: v.copy() for k, v in self.fields.items()}
