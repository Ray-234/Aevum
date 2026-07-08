# Climate Mechanism Modeling Plan

Status: active modeling plan
Last updated: 2026-07-07

This document defines the mechanism-first modeling work that must happen before
the next accepted climate implementation change.  It does not replace the
real-Earth single-subgraph replay route in
`docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md`.  The replay order remains:

1. R0 replay harness and invariants.
2. R1 boundary conditions and forcing.
3. R2a atmospheric pressure/source geometry.
4. R2b pressure-to-surface wind translation.
5. R3 ocean dynamics and basin coupling.
6. R4 SST and energy closure.
7. R5 moisture source and transport.
8. R6 seasonal precipitation and hydroclimate objects.
9. R7 cryosphere, cloud, and vegetation feedback.
10. R8 Koppen/classes and biomes.
11. R9 generated-world regression and promotion.

The current active replay packet remains R2a seasonal SLP / pressure-source
geometry.  Mechanism modeling is a design-and-contract step inside that route,
not a new fitting target and not a permission to tune downstream maps.

## Non-Negotiable Rules

- Fit one real-Earth subgraph at a time.  Generated worlds are inactive until
  R9 guardrail checks.
- Read maps before metrics.  Metrics can reject a bad run, but they cannot
  accept a geographically wrong map.
- Do not tune downstream products to hide upstream failures.  Wind, currents,
  SST, moisture, precipitation, ice, Koppen, and biomes remain blocked if R2a
  pressure/source geometry is visibly wrong.
- Parameter sweeps may only refine a named mechanism after the mechanism's
  inputs, outputs, budget checks, and map residual owner are explicit.
- The climate model should remain a reduced model, not a full GCM.  The target
  is interpretable physical causality: energy, momentum, water, and object
  layers that make Earth replay maps readable and extend to generated worlds.

## Hard M0/M1 Field Contract

Status, 2026-07-07:

- M0 uses existing shared geography primitives rather than a duplicate `geo.*`
  namespace.
- M0 now has a real-Earth replay override for major-ocean semantic basins:
  Atlantic, Pacific, Indian, Arctic, and Southern Ocean.  This is an Earth
  calibration support layer, not a generated-world shortcut.
- M1 now has diagnostic-only archive fields.  They are generated from existing
  state and do not change temperature, pressure, wind, ocean, moisture, or
  precipitation physics.
- All fields below are written to `terminal_climate_arrays.npz` when present.

| Layer | Field | Shape | Unit | Producer | Current source | Consumer | Mask/Range rule |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M0 | `terrain.elevation_m` | `(n,)` | m | terrain/reference replay | frozen Earth or generated terrain | M1, M2, M3, M6 | finite; sea level is 0 m in Earth replay |
| M0 | `climate.coast_distance` | `(n,)` | 1 | climate geography primitives | graph distance from coast | M1, M2, M5, M6 | normalized 0-1 |
| M0 | `climate.continent_id` | `(n,)` | id | climate geography primitives | land connected components | M2/M6 objects | ocean is -1 |
| M0 | `climate.continent_interiority` | `(n,)` | 1 | climate geography primitives | scaled land distance and continent size | M1, M2, monsoon support | land 0-1, ocean 0 |
| M0 | `ocean.basin_id` | `(n,)` | id | climate geography primitives | terrain basin ids, ocean connected components, or real-Earth major-ocean semantic ids | M2, M4, M5, M6 | ocean >=0, land -1 |
| M0 | `ocean.shelf_index` | `(n,)` | 1 | climate geography primitives | coast distance, depth, terrain shelf data | M4/M5 support | ocean 0-1, land 0 |
| M0 | `ocean.strait_index` | `(n,)` | 1 | climate geography primitives | narrow ocean adjacency | M4 gateway exchange | ocean 0-1, land 0 |
| M0 | `terrain.barrier_index` | `(n,)` | 1 | climate geography primitives | smoothed topographic height + gradient | M2, M3, M6 | land 0-1, ocean 0 |
| M0 | `terrain.wind_gap_index` | `(n,)` | 1 | climate geography primitives | low passes through barrier belts | M3/M6 passability | land 0-1, ocean 0 |
| M1 | `climate.seasonal_insolation_anomaly` | `(4,n)` | W m^-2 | climate energy diagnostic | seasonal TOA flux minus cell four-season mean | M1 map read, later pressure genesis | finite; seasonal mean per cell is 0 |
| M1 | `climate.surface_heat_capacity_class` | `(n,)` | 1 | climate energy diagnostic | ocean mask + continentality | M1/M2 attribution | ocean near 1, maritime land > interior land |
| M1 | `climate.land_thermal_anomaly` | `(4,n)` | K | climate energy diagnostic | seasonal land temperature minus cell annual mean | R2a land thermal low/high attribution | land finite, ocean 0 |
| M1 | `climate.ocean_mixed_layer_thermal_anomaly` | `(4,n)` | K | climate energy diagnostic | seasonal SST minus cell annual SST mean | R2a ocean pressure/SST attribution | ocean finite, land 0 |
| M1 | `climate.elevation_lapse_cooling` | `(n,)` | K | climate energy diagnostic | lapse rate times climate elevation | R2a/R6 terrain attribution | non-negative, clipped diagnostic |
| M1 | `climate.snow_ice_albedo_support` | `(4,n)` | 1 | climate energy diagnostic | seasonal sea ice over ocean, seasonal snow over land | R7 later; R2a polar residual attribution | 0-1 |
| M1 | `climate.sst_gradient_support` | `(4,n)` | 1 | climate energy diagnostic | normalized graph gradient of seasonal SST | R2a stationary waves, R3/R4 later | ocean 0-1.5, land 0 |
| M1 | `climate.same_latitude_sst_anomaly` | `(4,n)` | K | climate energy diagnostic | ocean SST minus same-latitude ocean band mean | over-zonal SST/pressure attribution | ocean finite, land 0 |
| M1 | `climate.land_sea_thermal_contrast` | `(4,n)` | K | climate energy diagnostic | locally spread land thermal anomaly minus ocean mixed-layer anomaly | R2a pressure-source attribution | finite, clipped to +/-45 K |

Current rendered evidence:

- `replay_m0_boundary_support_contact_sheet.png`
- `replay_m1_energy_support_contact_sheet.png`

Current R2a major-ocean checkpoint:

- `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/replay_m0_boundary_support_contact_sheet.png`
- Atlantic, Pacific, Indian, Arctic, and Southern Ocean ids are readable in the
  M0 support map.
- The M0 basin-id blocker is resolved.

Current M2 pressure-genesis checkpoint:

- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/real_earth_pressure_replay_contact_sheet.png`
- M2 v1 derives pressure-source objects from ocean-basin ids, open-ocean
  support, SST-front support, same-latitude SST anomaly, continent thermal
  anomaly, and terrain barrier support.
- M2 v2 separates causal source support from result diagnostics with signed
  source and per-source support maps, so map reading can distinguish generated
  Aleutian/Icelandic/Southern Ocean sources from inherited broad pressure
  support.
- M2 v4 adds Southern Ocean sector/wavenumber-front gating and continent-level
  thermal-center objects.  This is still a reduced pressure-source model, but
  causal source maps now show Aleutian/Icelandic source patches, Southern
  Ocean sectors, and continent-centered thermal refinements.
- M2 v6 adds a bounded source-to-pressure wave-transfer field,
  `atmosphere.pressure_genesis_wave_transfer`, so causal source objects can
  project into the final pressure proxy without changing R2b wind, R3 currents,
  SST, moisture, precipitation, ice, Koppen, biomes, or generated worlds.  The
  transfer is gated off over polar ice-cap edge geometry so Antarctic margin
  artifacts do not become ordinary continental stationary waves.
- M2 v7 replaces the v6 mostly isotropic transfer with a local weighted
  waveguide diffusion.  Ocean transfer follows SST-front support,
  same-latitude SST anomaly, open-ocean exposure, and basin/subpolar support;
  land transfer follows coast-strength, barrier, land-source, and terrain-wave
  support with the same non-polar gate.
- M2 v8 adds bounded object-level projection on top of the v7 waveguide:
  subpolar ocean-low support, land thermal-center source support, and
  terrain-wave source support are projected into the final pressure proxy with
  constrained footprints.  The transfer amplitude is still bounded and
  diagnostic-only for R2a; no downstream wind/current/SST/moisture/biome
  fitting is used.
- M2 v9 replaces the remaining ordinary weighted diffusion in the M2
  source-to-pressure transfer with directional diffusion.  Ocean transfer uses
  a storm-track/SST-front axis; land transfer uses terrain-barrier and coastal
  axes.  This implements the anisotropic-footprint contract without tuning
  downstream wind/current/SST/precipitation/biome fields.
- M2 v10 adds a bounded final-pressure expression solve.  It uses
  source-supported nonzonal pressure anomaly, aligned to the signed M2 source
  and transfer field, to let pressure-center objects affect the final pressure
  proxy more clearly without painting a whole latitude band.
- M2 v16 keeps the v10 final-expression solve and adds downwind cold-continent
  support to source generation.  Winter continental cold-high support is spread
  eastward/downwind into nearby open ocean and used to shape subpolar ocean-low
  source selection, so Aleutian/Icelandic-type lows are partly derived from
  cold continent plus westerly storm-track geometry rather than only from broad
  open-ocean latitude/front scoring.  Southern Hemisphere downwind support is
  damped so the Southern Ocean front/wavenumber gate remains primary.
- M2 v17 keeps the v16 source geometry but changes the final pressure
  synthesis.  Offline reconstruction of the v16 arrays showed that raw source
  strengthening degrades the SLP pattern, while stronger source-to-pressure
  transfer improves the Earth replay.  v17 therefore reduces the source's
  direct pressure contribution and increases the bounded transfer contribution;
  `atmosphere.pressure_genesis_source` stays a causal source diagnostic, while
  `atmosphere.pressure_genesis_wave_transfer` carries more of the actual
  pressure expression.
- M2 v18 keeps the v17 synthesis balance and adds a bounded transfer-morphology
  adjustment.  It deepens already source-supported ocean-low transfer cores and
  damps broad subtropical-high, land-core, and terrain-wave transfer
  contributions that made the final pressure proxy too smooth or fragmented.
  This is still M2-only and uses existing support fields; it does not add
  coordinate-specific targets or seasonal gain tables.
- M2 v19 keeps the v18 source/transfer/morphology structure and adds a bounded
  thermal-phase adjustment to the transfer expression.  Same-latitude SST and
  land-temperature anomalies decide whether high-support, low-support, and
  land-source footprints should strengthen or weaken the final pressure proxy,
  so cold ocean/land phases can support highs while warm phases do not promote
  the same support into the wrong seasonal sign.  This remains M2-only and does
  not tune downstream wind, currents, SST, precipitation, biomes, or generated
  worlds.
- M2 v25 keeps the v19 thermal-phase transfer, tightens the Southern Ocean
  wavenumber/front gate, and changes final pressure synthesis from a uniform
  direct-source weight to a domain-weighted source expression.  Land and
  Southern Ocean sources get less direct overprint; North Hemisphere mid/high
  latitude ocean-low cores retain more direct expression so Aleutian/Icelandic
  lows are not erased while the Southern Ocean band is softened.
- M2 v30 adds two source-placement refinements on top of v25: a post-normalized
  Northern Hemisphere winter coastal/front supplement for basin-edge subpolar
  lows, and a bounded shoulder-season warm-ocean/front low-source candidate.
  The shoulder candidate is component-selected but amplitude-limited after
  selection, so SON North Pacific/North Atlantic source objects can appear
  without turning MAM/SON into full winter-strength low bands.
- M2 v38 keeps the v30 shoulder-season source but makes it basin-scale and
  open-ocean weighted so SON North Pacific remains present while the
  North Atlantic shoulder low no longer over-deepens as strongly.  It also
  deepens Northern Hemisphere winter subpolar ocean-low expression and replaces
  the earlier Southern Ocean fixed-phase source emphasis with a front/shelf/SST
  support gate plus a signed Southern Ocean wave-transfer anomaly.
- M2 v40 keeps v38 and strengthens only the Northern Hemisphere winter
  subpolar ocean-low compact-expression term.  It deepens existing
  Aleutian/Icelandic-type cores without creating new MAM/SON/JJA low-source
  activation and improves ocean pressure correlation.
- M2 v44 keeps v40's compact NH winter low expression, adds a narrow coastal
  land inheritance term for DJF subpolar ocean lows, and adds a bounded
  poleward subpolar SST-front / storm-track low-support floor.  v43's attempt
  to put the lee-low term into object selection was rejected because it stole
  support from the Icelandic low and degraded ocean pressure correlation.
- M2 v45 keeps v44 and adds an Atlantic-Arctic gateway low-support floor.  It
  uses graph distance from existing Atlantic subpolar low support through the
  Atlantic/Arctic ocean domain, then only affects nearby 62-80 N Arctic-basin
  marginal seas with shelf and SST-front support.  This improves Nordic/
  Greenland/Barents winter pressure without applying a blanket Arctic low.
- M2 v46 keeps v45 and adds a Southern Ocean shoulder-season low-support floor.
  It uses the existing Southern Ocean front/shelf/wave gate to create MAM/SON
  low-source support in the semantic Southern Ocean basin, instead of leaving
  shoulder seasons with zero M2 source.
- R2a pressure is still not accepted.  The current source and transfer maps are
  auditable.  v46 is the current worktree checkpoint: it keeps the previously
  missing SON North Pacific source, improves DJF Aleutian/Icelandic expression,
  improves North Atlantic edge/Labrador and Nordic/Arctic gateway winter
  pressure, gives the Southern Ocean an explicit signed wave-transfer pattern,
  and adds missing Southern Ocean shoulder-season source support.  The final
  replay remains too smooth and not yet center-organized enough for R2a
  promotion.  The next owner remains M2 pressure genesis, specifically
  Southern Ocean wave-sector sharpening and residual North Atlantic compactness,
  not downstream tuning.

These are produced by `real-earth-wind-replay` when the replay arrays include
the fields above.  They are upstream evidence sheets for R2a; they are not
acceptance maps for R2a by themselves.

## Physical Coupling Basis

The climate system must be represented as a shared energy, momentum, and water
budget:

1. Boundary geometry:
   latitude, season, elevation, slope/aspect, land-sea mask, coast distance,
   ocean basin ids, gateways/straits, shelves/slopes, bathymetry, roughness,
   mountain barriers, and wind gaps.
2. Energy state:
   insolation, albedo, land/ocean heat capacity, elevation lapse-rate cooling,
   mixed-layer heat storage, SST gradients, sea-ice/snow support, and
   greenhouse background.
3. Pressure and wind:
   thermal lows/highs, subtropical highs, subpolar lows, ITCZ migration,
   stationary waves from land-sea contrast and mountains, pressure gradients,
   Coriolis turning, boundary-layer drag, terrain blocking, gap/channel winds,
   monsoon reversals, and storm-track jets.
4. Ocean circulation:
   wind stress, Ekman transport, wind-stress curl, basin gyres, western
   boundary currents, eastern boundary upwelling, equatorial upwelling,
   strait/gateway exchange, sea-ice edge effects, and reduced density
   overturning where temperature/salinity proxies exist.
5. SST and heat closure:
   shortwave input, ocean heat transport, upwelling/cold-tongue cooling,
   evaporation/latent heat loss, sensible heat exchange, sea-ice damping, and
   ocean-to-atmosphere feedback.
6. Moisture and precipitation:
   ocean evaporation, land evapotranspiration proxy, wind-routed moisture,
   convergence/divergence, ITCZ rainfall, storm-track precipitation, monsoon
   inflow, orographic windward rainfall, leeward rain shadow, convection over
   warm/moist land, snow processes, runoff, and soil/vegetation feedback.

## Submodel Contracts

### M0. Shared Boundary State

Purpose:

- Build a single reusable boundary-state bundle consumed by pressure, wind,
  ocean, moisture, and precipitation.

Required fields:

- `geo.latitude`, `geo.longitude`, `geo.cell_area`
- `geo.land_mask`, `geo.ocean_mask`, `geo.elevation_m`, `geo.bathymetry_m`
- `geo.slope`, `geo.aspect`, `geo.roughness`, `geo.barrier_strength`
- `geo.wind_gap_support`, `geo.coast_distance`, `geo.shelf_slope_class`
- `geo.ocean_basin_id`, `geo.gateway_id`, `geo.strait_support`
- `geo.continent_id`, `geo.catchment_id` when available

Map-read checks:

- Andes, Rockies, Himalaya/Tibet, East African Highlands, Greenland,
  Antarctica, major coastlines, ocean basins, and major straits must be visible
  in the relevant support maps.
- In real-Earth replay, Atlantic, Pacific, Indian, Arctic, and Southern Ocean
  support must be distinguishable before M2 ocean pressure centers are judged.
- Boundary maps are diagnostics only; they do not replace R2a pressure fitting.

### M1. Energy Boundary Layer

Purpose:

- Provide physically meaningful thermal drivers for pressure without tuning
  pressure directly.

Required fields:

- seasonal insolation anomaly
- surface heat-capacity class
- seasonal land thermal anomaly
- seasonal ocean mixed-layer thermal anomaly
- elevation lapse-rate cooling
- provisional snow/ice albedo support
- SST gradient support and same-latitude SST anomaly
- land-sea thermal contrast by coast/basin/continent object

Map-read checks:

- Land seasonal thermal response must be stronger than ocean response.
- High terrain must be colder than same-latitude low terrain.
- Ocean thermal anomalies should not be purely latitudinal where major currents
  or upwelling zones are present.

### M2. R2a Pressure Genesis

Purpose:

- Derive seasonal pressure/source geometry from energy state and boundary
  geometry, not from isolated painted blobs.

Required fields:

- `atmosphere.thermal_low_support`
- `atmosphere.thermal_high_support`
- `atmosphere.subtropical_high_support`
- `atmosphere.subpolar_low_support`
- `atmosphere.itcz_convergence_support`
- `atmosphere.land_sea_pressure_proxy`
- `atmosphere.stationary_wave_pressure_support`
- `atmosphere.pressure_center_support`

Required objects:

- pressure centers with season, sign, land/ocean/mixed domain, dominant
  continent/basin, centroid, area, strength, and source attribution.
- stationary-wave patches with source attribution: mountain, coast, SST front,
  land-sea thermal contrast, or basin geometry.

Map-read checks against real Earth:

- Aleutian/Icelandic winter lows should appear as ocean-basin-scale structures,
  not as generic latitude bands.
- Siberian/North American winter highs and Asian summer thermal lows should be
  tied to continental thermal state and terrain.
- Southern Ocean pressure should not collapse into a single broad zonal stripe.
- Coast and mountain stationary waves should be visible but not noisy.

## R2a Residual Owner Matrix

Use this matrix before any R2a pressure code change.  The owner named here
determines what may be changed.

| Visible R2a residual | First maps to inspect | Likely owner | Allowed action |
| --- | --- | --- | --- |
| Continental winter highs are missing or too weak | `land_thermal_anomaly`, `elevation_lapse_cooling`, `surface_heat_capacity_class`, Earth SLP | M1 if cold anomaly is absent; M2 if M1 is present but pressure ignores it | Fix M1 thermal driver only if M1 map is wrong; otherwise change pressure genesis weighting/object construction |
| Continental summer lows are broad blocky blobs | `land_thermal_anomaly`, `land_sea_thermal_contrast`, `continent_interiority`, pressure residual | M2 pressure genesis / smoothing | Replace blob-like response with continent/coast/terrain-supported pressure-center objects |
| Aleutian or Icelandic winter lows are missing | `same_latitude_sst_anomaly`, `sst_gradient_support`, `ocean.basin_id`, Earth SLP | M1 if SST/front support absent; M2 if support exists but no ocean pressure object | Do not tune wind; repair ocean-basin pressure genesis after confirming M1 support |
| Southern Ocean collapses into a zonal band | `ocean.basin_id`, `same_latitude_sst_anomaly`, `sst_gradient_support`, pressure zonal residual | M0/M1 if support is purely zonal; M2 if support is structured but ignored | Add basin/front/stationary-wave source geometry; do not compensate in R2b wind |
| Mountain/coast stationary waves are absent | `terrain.barrier_index`, `climate.coast_distance`, `climate.coast_strength`, `stationary_wave_pressure_support` | M0 if barrier/coast support is wrong; M2 if support exists but pressure object missing | Repair geography primitive or stationary-wave pressure support, not precipitation |
| Pressure centers exist but are too fragmented/noisy | `pressure_center_support`, `pressure_center_id`, M0/M1 support sheets | M2 object extraction | Adjust component/object extraction thresholds and continuity, not upstream energy |
| Pressure centers are smooth but displaced by one basin/continent | M0 basin/continent ids, M1 thermal/SST supports, pressure residual | M0 ownership if ids wrong; M2 ownership if ids right | Repair basin/continent primitive or pressure-center source attribution |
| R2b wind looks wrong but R2a pressure maps are wrong | R2a pressure contact sheet first | R2a upstream blocker | Do not tune R2b until R2a map-read acceptance |

### M3. R2b Surface Wind Translation

Purpose:

- Convert pressure/source geometry into near-surface winds with momentum
  constraints.

Required fields:

- pressure-gradient vectors
- Coriolis/geostrophic tendency
- boundary-layer drag and terrain roughness drag
- terrain blocking/gap-channel tendency
- divergence/convergence
- wind stress vector and wind-stress curl
- seasonal u10/v10 replay fields

Map-read checks against real Earth:

- Trades, westerlies, polar easterlies, monsoon reversals, and storm-track
  corridors must be readable as vector structures.
- Wind should bend around mountains and intensify through plausible gaps
  instead of crossing terrain as a latitude-only field.

### M4. R3 Ocean Dynamics

Purpose:

- Derive currents from wind stress and basin geometry.

Required fields:

- wind-stress curl
- basin streamfunction
- gyre id and gyre strength
- western/eastern boundary current support
- coastal/equatorial upwelling support
- downwelling support
- gateway/strait exchange
- reduced overturning support where temperature/salinity proxies are available

Map-read checks against real Earth:

- North Atlantic, North Pacific, South Atlantic, South Pacific, and Indian Ocean
  gyres must be basin-confined.
- Gulf Stream/Kuroshio-like western boundary currents should be narrow and
  strong relative to eastern boundary currents.
- Peru/Humboldt, Benguela, Canary/California-like eastern boundary upwelling
  zones should cool SST and affect moisture/biome later.

### M5. R4 SST and Heat Closure

Purpose:

- Solve SST as an energy balance rather than a post-hoc color field.

Required fields:

- net surface heat tendency
- current heat transport tendency
- upwelling cooling tendency
- evaporation/latent heat loss
- sea-ice damping/support
- final seasonal SST
- SST-front support

Map-read checks against real Earth:

- Western boundary warm tongues and eastern boundary cool tongues should be
  visible.
- Equatorial cold tongue and Southern Ocean gradients should be plausible.
- Same-latitude residual maps must be used to detect over-zonal SST.

### M6. R5-R6 Moisture and Precipitation

Purpose:

- Generate precipitation from moisture source, transport, convergence, and lift.

Required fields:

- evaporation source strength by basin
- moisture flux vector
- moisture convergence/divergence
- source-basin attribution for receiver cells
- orographic uplift and rain-shadow support
- ITCZ/storm-track/monsoon precipitation support
- seasonal precipitation

Required objects:

- ITCZ rain belts
- monsoon inflow corridors
- storm-track precipitation corridors
- orographic windward rain regions
- rain-shadow regions
- receiver catchments with source-basin ledgers

Map-read checks against real Earth:

- Amazon/Congo/Maritime Continent wet regions must come from tropical moisture
  convergence, not arbitrary biome thresholds.
- Sahara/Arabia/Australia dry belts must be tied to subtropical subsidence and
  moisture-path limits.
- Andes, Cascades/Sierra Nevada, Himalaya/Tibet, New Zealand, and island chains
  must show windward/leeward contrast.

## Earth Replay Evidence Packet

Every implementation packet must write:

- real-Earth reference map for the active subgraph
- Aevum replay map on the same grid
- residual or vector-error map
- relevant support maps from upstream fields
- object overlay/contact sheet when the active layer has process objects
- map-read attribution note naming real structures, replay successes, replay
  failures, and residual owner
- numeric regression table after map reading, not before

## Microbenchmarks

M0 boundary state:

- field presence, finite values, mask consistency, area-weighted coverage
- coast/basin/gateway/barrier object continuity

M1 energy:

- land-ocean seasonal amplitude contrast
- elevation lapse-rate monotonic sanity by same-latitude bins
- same-latitude SST anomaly non-zonality guard

M2 pressure:

- pressure-center object count and domain balance by season
- support overlap with real Earth pressure-center basins
- excessive zonal-band residual guard
- stationary-wave support continuity near major terrain/coast triggers

M3 wind:

- vector angular error map
- divergence/convergence placement
- terrain-crossing penalty
- monsoon reversal and storm-track corridor placement

M4 ocean:

- basin-confined streamfunction continuity
- western-boundary-current narrowness/strength
- upwelling placement along eastern boundaries and equator
- no-current-through-land invariant

M5 SST:

- same-latitude SST residual
- western/eastern boundary contrast
- upwelling cold-tongue placement
- bounded heat-budget residual

M6 moisture/precipitation:

- source-to-receiver moisture ledger closure
- wet/dry region support attribution
- windward/leeward precipitation contrast
- monsoon wet-season timing

## Reference Case Studies

Use these real-Earth structures as map-reading anchors:

- South Asian monsoon: land-sea thermal contrast, Tibetan Plateau influence,
  Indian Ocean moisture source, summer inflow, orographic precipitation.
- North American monsoon: summer continental heating, Gulf of California/Gulf of
  Mexico moisture pathways, terrain channeling.
- Andes / Humboldt / Atacama: mountain rain shadow, cold eastern boundary
  current, coastal upwelling, subtropical subsidence.
- Gulf Stream / North Atlantic storm track: western boundary current, SST front,
  storm-track precipitation, Icelandic Low.
- Kuroshio / North Pacific storm track: western boundary current, SST front,
  Aleutian Low.
- Southern Ocean: circumpolar westerlies, open-ocean storm belt, upwelling, sea
  ice edge.
- Maritime Continent: complex land-sea mask, warm pool evaporation, convection,
  seasonal monsoon modulation.

## Next Implementation Entry

Current status as of the R2a M2 v63 checkpoint:

- The M0/M1 diagnostic contract exists and is active for real-Earth pressure
  replay.  Each R2a iteration now renders Earth SLP, Aevum pressure, residual,
  M0 support, M1 support, M2 source, and M2 source-to-pressure transfer maps.
- Major-ocean semantic basin support exists for Atlantic, Pacific, Indian,
  Arctic, and Southern Ocean in Earth replay.
- M2 now contains object/source support for Aleutian/Icelandic winter lows,
  Atlantic-Arctic gateway winter lows, continental seasonal pressure centers,
  Northern Hemisphere MAM/JJA subpolar-ocean highs, a MAM Arctic freeze-ocean
  high object, JJA North Pacific / Gulf Alaska / central Pacific high objects,
  a JJA Eurasian summer thermal-low belt, a JJA western Pacific marginal-sea
  trough, land shoulder-season phase correction, MAM North America spring and
  plains land-high objects, Canadian / central-west Arctic freeze-high objects,
  DJF North America winter-high relief, SON boreal autumn land-high relief,
  SON North Atlantic autumn low support, and Southern Ocean shoulder-season
  source plus latitude / sector wave-transfer terms.
- The current evidence packet is v63:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/`.
- v63 improves standardized-pressure MAE all/land/ocean from v52
  `0.284/0.394/0.239` to `0.271/0.373/0.229`, all/land/ocean pressure
  correlation from `0.626/0.689/0.495` to `0.649/0.699/0.551`, MAM
  standardized-pressure MAE / zonal-anomaly correlation from `0.313/0.566`
  to `0.290/0.615`, JJA standardized-pressure MAE / zonal-anomaly correlation
  from `0.207/0.734` to `0.192/0.778`, and SON standardized-pressure MAE from
  `0.388` at v60 to `0.378`.
- R2a is still not accepted.  v52-v63 fix the largest land shoulder-season,
  MAM Arctic/Greenland/Beaufort ocean, JJA North Pacific/Gulf Alaska, DJF
  Atlantic-Arctic / Aleutian winter-low, DJF North America winter over-high,
  SON boreal autumn land over-high, SON North Atlantic / Icelandic autumn
  under-low, the v61 Southern Ocean side effect, and the JJA western Pacific
  sign error.  Remaining visible owners are residual MAM Baffin-Labrador /
  central Arctic under-high, JJA North Atlantic / Icelandic under-high,
  residual JJA Japan / East China Sea and Arabia / Iran summer low
  underexpression, plus residual SON high-latitude texture / Antarctica edge
  artifacts.

The next code step should still be R2a/M2-only:

1. Read the v63 Earth SLP, replay pressure, standardized residual, zonal
   residual, source, transfer, M0, and M1 maps before changing code.
2. Select the next owner by map read: residual MAM Baffin-Labrador / central
   Arctic under-high, JJA North Atlantic / Icelandic under-high, JJA Japan /
   East China Sea and Arabia / Iran summer low underexpression, or SON
   high-latitude texture.
3. Preserve the v49 Southern Ocean latitude split, v48 MAM latitude narrowing,
   v50 Atlantic-Arctic gateway transfer, v51 Beaufort-protected JJA subpolar
   high support, v52 land shoulder-season phase correction, v53 MAM Arctic
   freeze-ocean high support, v54/v60 JJA North Pacific high support, v55 MAM
   North America spring land-high support, v58 DJF Atlantic-Arctic / Aleutian
   low plus North America winter-high relief improvements, v59 MAM
   central/eastern North America plains-high plus Canadian Arctic support, and
   v62 SON autumn land / North Atlantic / Southern Ocean corrections, and v63
   JJA Eurasian thermal-low / western Pacific trough corrections.
4. Keep R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated
   worlds blocked as fitting targets until R2a passes map-read acceptance.

## Sources

- NOAA National Ocean Service, trade winds and subtropical highs:
  https://oceanservice.noaa.gov/facts/tradewinds.html
- NOAA/NESDIS SciJinks, monsoon land-ocean seasonal thermal contrast:
  https://scijinks.gov/what-is-a-monsoon/
- NOAA Climate.gov, annual migration of the tropical rain belt:
  https://www.climate.gov/news-features/understanding-climate/annual-migration-tropical-rain-belt
- NOAA National Ocean Service, windward/leeward orographic effect:
  https://oceanservice.noaa.gov/facts/windward-leeward.html
- NOAA National Ocean Service, boundary currents and gyres:
  https://oceanservice.noaa.gov/education/tutorial_currents/04currents3.html
- NOAA National Ocean Service, Ekman spiral:
  https://oceanservice.noaa.gov/education/tutorial_currents/04currents4.html
- NOAA National Ocean Service, upwelling:
  https://oceanservice.noaa.gov/facts/upwelling.html
- NOAA Ocean Exploration, ocean currents:
  https://oceanexplorer.noaa.gov/ocean-fact/currents/
- NOAA Ocean Exploration, ocean effects on climate and weather:
  https://oceanexplorer.noaa.gov/ocean-fact/climate/
- NOAA/NESDIS SciJinks, jet stream and storm placement:
  https://www.nesdis.noaa.gov/about/k-12-education/atmosphere/what-the-jet-stream
- IPCC AR6 WGI Chapter 8, water cycle and circulation:
  https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-8/
- IPCC AR6 WGI Chapter 9, ocean circulation and heat transport:
  https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-9/
