"""PlanetSpec: the input parameter specification.

Differences between worlds must come from THREE levels (project plan section 5):

  1. parameter differences  -- this dataclass
  2. regime differences     -- e.g. tectonic mode (not all habitable worlds have
                               modern Earth-style mobile-lid plate tectonics)
  3. history differences    -- emergent from impacts, volcanism, rifting, biology

The engine's internal clock is "time since formation"; ``geologic_labels`` is only
a display mapping for the Earth preset.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional

from aevum.core.units import CONSTANTS


class TectonicRegime(str, Enum):
    """Tectonic modes are a continuum/discrete set, not a single mode."""

    MOBILE_LID = "mobile_lid"          # modern Earth-style plate tectonics
    SLUGGISH_LID = "sluggish_lid"      # slow, partial recycling
    EPISODIC_LID = "episodic_lid"      # intermittent overturn
    STAGNANT_LID = "stagnant_lid"      # volcanic-plain dominated, no plates


@dataclass
class StellarSpec:
    mass_solar: float = 1.0
    luminosity_solar_now: float = 1.0      # present luminosity (L_sun)
    spectral_temp_k: float = 5772.0
    age_at_start_gyr: float = 0.0          # star age when planet forms


@dataclass
class OrbitSpec:
    semi_major_axis_au: float = 1.0
    eccentricity: float = 0.0167
    obliquity_deg: float = 23.44
    rotation_period_hours: float = 24.0
    tidally_locked: bool = False


@dataclass
class CompositionSpec:
    mass_earth: float = 1.0
    radius_earth: float = 1.0
    core_mass_fraction: float = 0.32
    radiogenic_abundance: float = 1.0     # relative to Earth (U, Th, K)
    water_inventory_earth: float = 1.0    # total H2O budget vs Earth oceans
    initial_internal_temp_k: float = 3500.0


@dataclass
class AtmosphereSpec:
    surface_pressure_bar: float = 1.0
    co2_bar: float = 280e-6               # partial pressure of CO2 (bar)
    o2_fraction: float = 0.0              # starts anoxic by default
    n2_fraction: float = 0.78


@dataclass
class BiosphereSpec:
    life_origin_myr: Optional[float] = 500.0   # None = lifeless world
    allow_oxygenic_photosynthesis: bool = True
    allow_land_colonization: bool = True


@dataclass
class PlanetSpec:
    name: str = "unnamed"
    seed: int = 1
    t_end_myr: float = 4500.0             # internal age, not "Quaternary"

    stellar: StellarSpec = field(default_factory=StellarSpec)
    orbit: OrbitSpec = field(default_factory=OrbitSpec)
    composition: CompositionSpec = field(default_factory=CompositionSpec)
    atmosphere: AtmosphereSpec = field(default_factory=AtmosphereSpec)
    biosphere: BiosphereSpec = field(default_factory=BiosphereSpec)

    initial_tectonic_regime: TectonicRegime = TectonicRegime.MOBILE_LID
    crust_strength: float = 1.0           # relative lithospheric strength
    target_land_fraction: float = 0.29    # initial continental crust fraction

    n_plates: int = 12
    impact_flux_scale: float = 1.0        # relative to Earth's bombardment history
    n_moons: int = 1

    grid_cells: int = 20000

    notes: str = ""

    # -- derived quantities ---------------------------------------------
    @property
    def radius_m(self) -> float:
        return self.composition.radius_earth * CONSTANTS.EARTH_RADIUS

    @property
    def mass_kg(self) -> float:
        return self.composition.mass_earth * CONSTANTS.EARTH_MASS

    @property
    def surface_gravity(self) -> float:
        return CONSTANTS.G * self.mass_kg / (self.radius_m ** 2)

    @property
    def surface_area_m2(self) -> float:
        import math
        return 4.0 * math.pi * self.radius_m ** 2

    def to_dict(self) -> dict[str, Any]:
        from dataclasses import asdict
        d = asdict(self)
        d["initial_tectonic_regime"] = self.initial_tectonic_regime.value
        return d
