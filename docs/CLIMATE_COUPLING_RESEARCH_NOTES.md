# Climate Coupling Research Notes

Status: active mechanism notes
Last updated: 2026-07-07

These notes support the real-Earth calibration track.  They summarize the
physical coupling among terrain, winds, ocean currents, temperature, and
precipitation, and explain why the next climate work should be mechanism-first
rather than parameter tuning.

## Core Coupling Graph

The climate pipeline should be treated as a coupled energy, momentum, and water
system:

1. Boundary conditions:
   terrain/elevation, land-sea mask, basin geometry, coastlines, gateways,
   roughness/barriers, latitude, season, rotation, and radiative forcing.
2. Energy state:
   insolation, land/ocean heat capacity, elevation lapse-rate cooling, ocean
   heat storage, ocean heat transport, SST, snow/ice/albedo, and greenhouse
   forcing.
3. Pressure and wind:
   pressure gradients from uneven heating, Coriolis turning, boundary-layer
   drag, terrain blocking, channeling, mountain waves, monsoon reversals, and
   storm-track jets.
4. Ocean dynamics:
   wind stress, Ekman transport, gyres, boundary currents, upwelling,
   downwelling, strait exchange, and density-driven overturning from
   temperature/salinity.
5. Moisture and precipitation:
   evaporation from SST and wind, moisture transport by circulation, convergence
   in ITCZ/storm tracks/monsoons, orographic uplift, rain shadow, runoff, soil
   moisture, vegetation, snow, and ice feedback.

The practical implication for Aevum is that R2 pressure/wind, R3 currents, R4
SST/energy, and R5-R6 moisture/precipitation cannot be made Earthlike by
adjusting independent scalar maps.  Each layer must expose process objects and
budgets that the next layer consumes.

## Physical Principles

Atmospheric circulation:

- Uneven solar heating creates pressure gradients and convection; Earth's
  rotation turns large-scale flow through the Coriolis effect.
- Subtropical high-pressure belts near 30 degrees feed trades and westerlies.
  The ITCZ is a convergence/rising-motion band, not just a latitude stripe.
- Monsoons are seasonal circulation reversals caused by land-sea thermal
  contrast, modified by elevated landmasses, surrounding oceans, and moisture
  availability.

Terrain:

- Elevation cools air through lapse-rate effects and changes snow/ice/albedo.
- Mountains force moist air to rise, cool, condense, and precipitate on
  windward slopes; leeward sinking produces rain shadows.
- Orography also modifies the wind field through blocking, channeling, lee
  troughing, mountain waves, barrier jets, and gap winds.
- Model resolution matters: coarse topography weakens orographic precipitation,
  local wind corridors, and regional water-cycle structure.

Ocean currents:

- Winds drive the upper ocean through friction and Coriolis-deflected Ekman
  transport.
- Basin geometry, coastlines, rotation, and wind-stress curl create gyres and
  boundary currents.
- Temperature and salinity control density-driven overturning and deep currents.
- Ocean currents redistribute heat and moisture sources, affecting SST,
  evaporation, atmospheric pressure, storm tracks, and precipitation.

Temperature and moisture:

- SST and land temperature determine sensible/latent heat fluxes and therefore
  pressure, wind, evaporation, and boundary-layer humidity.
- Warmer air can hold more water vapour, but precipitation location is governed
  by circulation and convergence, not temperature alone.
- Land usually has lower heat capacity and less moisture availability than the
  ocean, so land-sea thermal contrast and relative-humidity changes matter for
  monsoons and continental dryness.

Precipitation:

- Precipitation needs a moisture source, transport path, lifting/convergence
  mechanism, and sufficient thermodynamic support.
- Major process families should be separate objects: ITCZ rain belt, storm-track
  frontal precipitation, monsoon inflow, orographic windward rain, rain shadow,
  convection over warm/moist land, and polar/snow processes.

## Implications for Current Code

The current R2a pressure proxy has useful diagnostics but remains too map-like:

- It has broad land thermal lows/highs and broad ocean pressure sources.
- It does not yet solve pressure as the result of a closed energy/momentum
  balance.
- Its ocean basin centers and stationary waves are diagnostic objects, not
  causal objects.

Therefore, the next work should not be another local pressure-source parameter
change.  It should first refactor the replay logic into explicit sub-models:

1. R1/R4 energy boundary layer:
   derive land and ocean thermal state, SST gradients, elevation cooling,
   snow/ice albedo, and surface heat capacity as reusable fields.
2. R2a pressure genesis:
   derive pressure centers from thermal gradients, SST gradients, orography,
   land-sea contrast, and existing latitude circulation support.
3. R2b wind translation:
   derive near-surface winds from pressure gradients plus geostrophic balance,
   boundary-layer drag, terrain blocking, and gap/channel effects.
4. R3 ocean circulation:
   derive wind-driven currents from wind stress, basin geometry, Coriolis,
   coastlines, gateways, and upwelling/downwelling; add reduced density
   overturning only after SST/salinity proxies exist.
5. R4 SST/heat closure:
   update SST from insolation, ocean heat transport, upwelling, evaporation,
   sea ice, and land/ocean heat exchange.
6. R5-R6 moisture and precipitation:
   solve evaporation, moisture transport, convergence/lift, orographic
   precipitation, rain shadows, monsoon regions, and storm-track corridors as
   process objects with local conservation checks.

## Source Pointers

- NOAA National Ocean Service, trade winds and subtropical highs:
  https://oceanservice.noaa.gov/facts/tradewinds.html
- NOAA National Ocean Service, horse latitudes / subtropical high-pressure dry
  belts:
  https://oceanservice.noaa.gov/facts/horse-latitudes.html
- NOAA National Ocean Service, Ekman spiral:
  https://oceanservice.noaa.gov/education/tutorial_currents/04currents4.html
- NOAA National Ocean Service, thermohaline circulation:
  https://oceanservice.noaa.gov/education/tutorial_currents/05conveyor1.html
- NOAA Ocean Exploration, ocean effects on climate and weather:
  https://oceanexplorer.noaa.gov/ocean-fact/climate/
- NOAA National Weather Service glossary, orographic precipitation:
  https://forecast.weather.gov/glossary.php?word=RA
- IPCC AR6 WGI Chapter 8, water-cycle changes:
  https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-8/
- IPCC AR6 WGI Chapter 9, ocean, cryosphere, and sea level:
  https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-9/

