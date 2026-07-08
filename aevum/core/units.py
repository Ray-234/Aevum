"""Units, physical dimensions, constants and conserved quantities.

This is intentionally lightweight: full dimensional algebra is overkill for the
engine.  What we need is (a) every feature carries a canonical unit string so
nothing is dimensionless-by-accident, and (b) an enumerated set of globally
conserved quantities that the validation layer can audit.
"""
from __future__ import annotations

from enum import Enum


class Dimension(str, Enum):
    """Canonical physical dimension tags used in the feature registry."""

    DIMENSIONLESS = "dimensionless"
    LENGTH = "length"
    AREA = "area"
    VOLUME = "volume"
    TIME = "time"
    MASS = "mass"
    TEMPERATURE = "temperature"
    PRESSURE = "pressure"
    ENERGY = "energy"
    POWER = "power"
    ENERGY_FLUX = "energy_flux"          # W / m^2
    VELOCITY = "velocity"
    MASS_FLUX = "mass_flux"              # kg / m^2 / s
    CONCENTRATION = "concentration"
    MOLES = "moles"
    FRACTION = "fraction"
    ANGLE = "angle"


#: Canonical SI-ish unit string per dimension (informational / for display).
CANONICAL_UNIT: dict[Dimension, str] = {
    Dimension.DIMENSIONLESS: "1",
    Dimension.LENGTH: "m",
    Dimension.AREA: "m^2",
    Dimension.VOLUME: "m^3",
    Dimension.TIME: "s",
    Dimension.MASS: "kg",
    Dimension.TEMPERATURE: "K",
    Dimension.PRESSURE: "Pa",
    Dimension.ENERGY: "J",
    Dimension.POWER: "W",
    Dimension.ENERGY_FLUX: "W m^-2",
    Dimension.VELOCITY: "m s^-1",
    Dimension.MASS_FLUX: "kg m^-2 s^-1",
    Dimension.CONCENTRATION: "mol m^-3",
    Dimension.MOLES: "mol",
    Dimension.FRACTION: "1",
    Dimension.ANGLE: "rad",
}


class Conserved(str, Enum):
    """Globally conserved quantities audited by the validation layer.

    Conservation is checked as: total inventory change over a step must equal
    the net of declared sources minus sinks (within tolerance).
    """

    MASS = "mass"
    ENERGY = "energy"
    WATER = "water"
    SALT = "salt"
    CARBON = "carbon"
    OXYGEN = "oxygen"
    NITROGEN = "nitrogen"
    PHOSPHORUS = "phosphorus"
    SULFUR = "sulfur"
    CRUST_AREA = "crust_area"


class CONSTANTS:
    """Physical constants (SI)."""

    G = 6.67430e-11                 # gravitational constant
    STEFAN_BOLTZMANN = 5.670374419e-8
    SOLAR_LUMINOSITY = 3.828e26     # W
    SOLAR_MASS = 1.98892e30         # kg
    SOLAR_RADIUS = 6.957e8          # m
    AU = 1.495978707e11             # m
    EARTH_MASS = 5.972e24           # kg
    EARTH_RADIUS = 6.371e6          # m
    EARTH_GRAVITY = 9.80665         # m / s^2
    SECONDS_PER_YEAR = 31557600.0   # Julian year
    SECONDS_PER_MYR = 31557600.0 * 1e6
    R_GAS = 8.314462618             # J / mol / K
    LATENT_HEAT_VAPOR = 2.5e6       # J / kg
    LATENT_HEAT_FUSION = 3.34e5     # J / kg
    WATER_DENSITY = 1000.0          # kg / m^3
    ICE_DENSITY = 917.0
    CRUST_DENSITY_CONT = 2700.0     # kg / m^3
    CRUST_DENSITY_OCEAN = 2900.0
    MANTLE_DENSITY = 3300.0
    ZERO_C = 273.15                 # K


def myr_to_seconds(myr: float) -> float:
    return myr * CONSTANTS.SECONDS_PER_MYR
