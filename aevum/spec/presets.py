"""Six baseline worlds (project plan, Phase 0).

These exist so the engine is validated against a *combination* of worlds -- Earth
is only one test scenario and must never become the sole calibration target.
"""
from __future__ import annotations

from aevum.spec.planet_spec import (
    AtmosphereSpec,
    BiosphereSpec,
    CompositionSpec,
    OrbitSpec,
    PlanetSpec,
    StellarSpec,
    TectonicRegime,
)


def _earthlike() -> PlanetSpec:
    return PlanetSpec(
        name="earthlike_mobile_lid",
        seed=42,
        t_end_myr=4500.0,
        stellar=StellarSpec(mass_solar=1.0, luminosity_solar_now=1.0),
        orbit=OrbitSpec(semi_major_axis_au=1.0, eccentricity=0.0167,
                        obliquity_deg=23.44, rotation_period_hours=24.0),
        composition=CompositionSpec(water_inventory_earth=1.0),
        atmosphere=AtmosphereSpec(surface_pressure_bar=1.0, co2_bar=3e-3),
        biosphere=BiosphereSpec(life_origin_myr=500.0),
        initial_tectonic_regime=TectonicRegime.MOBILE_LID,
        target_land_fraction=0.29,
        n_plates=12,
        notes="Modern-Earth analogue; the reference (not the only) calibration world.",
    )


def _waterworld() -> PlanetSpec:
    return PlanetSpec(
        name="waterworld",
        seed=7,
        t_end_myr=4500.0,
        composition=CompositionSpec(water_inventory_earth=12.0),
        atmosphere=AtmosphereSpec(surface_pressure_bar=1.5, co2_bar=5e-3),
        biosphere=BiosphereSpec(life_origin_myr=400.0, allow_land_colonization=False),
        initial_tectonic_regime=TectonicRegime.MOBILE_LID,
        target_land_fraction=0.03,
        n_plates=10,
        notes="Deep global ocean, negligible exposed land.",
    )


def _arid() -> PlanetSpec:
    return PlanetSpec(
        name="arid_world",
        seed=101,
        t_end_myr=4500.0,
        composition=CompositionSpec(water_inventory_earth=0.15),
        atmosphere=AtmosphereSpec(surface_pressure_bar=0.6, co2_bar=8e-3),
        biosphere=BiosphereSpec(life_origin_myr=800.0),
        initial_tectonic_regime=TectonicRegime.MOBILE_LID,
        target_land_fraction=0.78,
        n_plates=9,
        notes="Water-poor; small seas, vast continents.",
    )


def _stagnant_lid() -> PlanetSpec:
    return PlanetSpec(
        name="stagnant_lid_world",
        seed=303,
        t_end_myr=4500.0,
        composition=CompositionSpec(radiogenic_abundance=0.7,
                                    initial_internal_temp_k=3200.0),
        atmosphere=AtmosphereSpec(surface_pressure_bar=2.0, co2_bar=0.1),
        biosphere=BiosphereSpec(life_origin_myr=1500.0,
                                allow_oxygenic_photosynthesis=True),
        initial_tectonic_regime=TectonicRegime.STAGNANT_LID,
        crust_strength=3.0,
        target_land_fraction=0.55,
        n_plates=1,
        notes="No plate recycling; volcanic plains dominate.",
    )


def _tidally_locked() -> PlanetSpec:
    return PlanetSpec(
        name="tidally_locked_world",
        seed=909,
        t_end_myr=4500.0,
        stellar=StellarSpec(mass_solar=0.4, luminosity_solar_now=0.04,
                            spectral_temp_k=3400.0),
        orbit=OrbitSpec(semi_major_axis_au=0.15, eccentricity=0.02,
                        obliquity_deg=2.0, rotation_period_hours=24 * 18,
                        tidally_locked=True),
        composition=CompositionSpec(water_inventory_earth=0.8),
        atmosphere=AtmosphereSpec(surface_pressure_bar=1.2, co2_bar=1e-2),
        biosphere=BiosphereSpec(life_origin_myr=900.0),
        initial_tectonic_regime=TectonicRegime.SLUGGISH_LID,
        target_land_fraction=0.4,
        n_plates=6,
        notes="M-dwarf, permanent day/night hemispheres.",
    )


def _snowball() -> PlanetSpec:
    return PlanetSpec(
        name="frozen_world",
        seed=555,
        t_end_myr=4500.0,
        stellar=StellarSpec(mass_solar=0.9, luminosity_solar_now=0.75),
        orbit=OrbitSpec(semi_major_axis_au=1.25, eccentricity=0.05,
                        obliquity_deg=18.0, rotation_period_hours=22.0),
        composition=CompositionSpec(water_inventory_earth=1.2),
        atmosphere=AtmosphereSpec(surface_pressure_bar=0.8, co2_bar=1e-3),
        biosphere=BiosphereSpec(life_origin_myr=1200.0),
        initial_tectonic_regime=TectonicRegime.MOBILE_LID,
        target_land_fraction=0.35,
        n_plates=11,
        notes="Cold orbit; prone to repeated global glaciation.",
    )


PRESETS = {
    "earthlike": _earthlike,
    "waterworld": _waterworld,
    "arid": _arid,
    "stagnant_lid": _stagnant_lid,
    "tidally_locked": _tidally_locked,
    "frozen": _snowball,
}


def get_preset(name: str) -> PlanetSpec:
    if name not in PRESETS:
        raise KeyError(f"unknown preset '{name}'. options: {sorted(PRESETS)}")
    return PRESETS[name]()
