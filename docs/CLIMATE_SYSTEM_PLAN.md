# Aevum Climate System Plan

Status: active planning document
Owner: climate / terrain / compiler integration
Last updated: 2026-07-07

This document archives the climate redesign plan and should be updated during
future climate work.  It is intentionally more detailed than the README roadmap:
the README describes the project; this file tracks the climate subsystem.

## Current Authoritative Workflow

The active development mode is real-Earth replay fitting, one subgraph at a
time.  For R0-R8, do not use generated virtual worlds to choose parameters or
judge mechanism changes.  The only acceptance evidence for an active phase is
the same-grid real-Earth reference subgraph, the same-grid Aevum real-Earth
replay subgraph, and their residual/error maps, followed by explicit
geoscience attribution.

Each iteration must:

1. Select one Earth reference subgraph, for example R2a seasonal SLP /
   pressure-source geometry.
2. Render the Earth reference map, replay map, and residual/vector-error map.
3. Read the geography first: name the real structures, replay successes, and
   replay failures.
4. Assign each major residual to an upstream owner or the active mechanism.
5. Change only that owner, then re-render the same Earth subgraph.

Global means, percentile envelopes, correlations, and gate tables are
regression tools.  They can reject a bad run, but they cannot promote a
subgraph whose map still has the wrong geography.  Downstream maps, including
currents, SST, precipitation, sea ice, Koppen classes, biomes, and generated
worlds, are observer-only until their upstream Earth replay layers pass.

## Mechanism-First Reset, 2026-07-07

The latest R2a pressure-source experiment showed that the next important task is
not parameter tuning.  A local basin-pressure-source tweak was tested in
`out_real_earth_pressure_replay_r2a_basin_pressure_v1_20260707/` and worsened
the ocean pressure correlation while producing no useful map-level repair; the
related code change was rejected and rolled back.  The partial v2 replay under
`out_real_earth_climate_replay_r2a_basin_pressure_v2_20260707/` is also a
rejected tuning artifact, not an accepted checkpoint.

Before another pressure, wind, current, SST, moisture, or precipitation code
change, the climate track must first define the coupled mechanism contracts:
terrain and land-sea geometry -> energy state -> pressure/wind -> ocean
circulation -> SST/heat closure -> moisture transport -> precipitation
objects.  The supporting research note is archived in
`docs/CLIMATE_COUPLING_RESEARCH_NOTES.md`, and the executable modeling plan is
archived in `docs/CLIMATE_MECHANISM_MODELING_PLAN.md`.

Immediate consequence:

- Keep R2a seasonal SLP / pressure-source geometry as the active replay target.
- Pause local pressure-source parameter work until the energy, momentum, and
  water-budget handoff fields are specified.
- Treat pressure-center and stationary-wave maps as diagnostics unless a later
  mechanism explicitly consumes them as causal objects.
- Continue using real-Earth subgraph maps as the fitting target; generated
  worlds remain R9 guardrails.

2026-07-07 implementation checkpoint:

- Added diagnostic-only M1 energy-boundary fields to `ClimateModule` and the
  terminal/replay archive: seasonal insolation anomaly, reduced surface
  heat-capacity class, land thermal anomaly, ocean mixed-layer thermal anomaly,
  elevation lapse cooling, snow/ice albedo support, SST-gradient support,
  same-latitude SST anomaly, and land-sea thermal contrast.
- Extended `real-earth-wind-replay` to render M0/M1 support contact sheets when
  these fields are present:
  `replay_m0_boundary_support_contact_sheet.png` and
  `replay_m1_energy_support_contact_sheet.png`.
- This is a diagnostic contract checkpoint only.  It does not change accepted
  pressure, wind, ocean, temperature, moisture, or precipitation mechanics.
- First visual attribution is archived in
  `docs/R2A_M0_M1_MAP_READ_ATTRIBUTION_20260707.md`.
- Added real-Earth major-ocean semantic basin support for R2a replay.  In
  Earth replay, `ocean.basin_id` now separates Atlantic, Pacific, Indian,
  Arctic, and Southern Ocean support instead of collapsing the ocean into one
  connected component.  Evidence:
  `out_real_earth_climate_replay_r2a_major_ocean_basins_20260707/` and
  `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/`.
- R2a remains unaccepted.  The M0 basin-id blocker is resolved, but M2
  pressure genesis still turns available energy/boundary support into
  over-smooth continental blobs and zonal ocean bands.
- Added M2 pressure-genesis v1: an object-based pressure-source refinement that
  extracts winter subpolar ocean-low candidates and weak subtropical highs from
  basin/front/open-ocean support, with small continent/terrain stationary-wave
  refinements.  It updates `atmosphere.land_sea_pressure_proxy`; it does not
  change R2b wind translation.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v1_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v1_20260707/`.
- M2 v1 is a checkpoint, not acceptance.  Map read improved North Pacific /
  North Atlantic low readability and introduced Southern Ocean segmentation,
  but pressure remains too smooth and partly over-zonal.
- Added M2 pressure-genesis v2: the causal pressure-source increment is now
  archived separately from the pressure-center diagnostic result as
  `atmosphere.pressure_genesis_source`,
  `atmosphere.ocean_pressure_low_source_support`,
  `atmosphere.ocean_pressure_high_source_support`,
  `atmosphere.land_pressure_source_support`, and
  `atmosphere.terrain_pressure_wave_source_support`.  Real-Earth pressure
  replay renders these maps in the pressure evidence packet.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v2_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/`.
- M2 v2 remains unaccepted but is cleaner than v1: causal sources are readable,
  ocean pressure correlation remains materially better than the M0 checkpoint,
  and MAE is close to the M0 checkpoint instead of the more degraded v1.
- Added M2 pressure-genesis v4: Southern Ocean low-source extraction now uses
  basin-sector labels plus a wavenumber/front gate, so the causal source reads
  as circumpolar pressure-wave sectors rather than one uniform annular band.
  Continental pressure refinements now use continent-level thermal-center
  objects and terrain/land-thermal-gradient support instead of another broad
  all-land scalar.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v4_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/`.
- M2 v4 is the current R2a checkpoint, not acceptance.  It improves ocean and
  land pressure correlations while keeping MAE near the M0 checkpoint, but the
  final pressure field is still smoother than real Earth SLP.

## Current State

The current climate module is a reduced seasonal model:

- Annual temperature is solved with a diffusive energy-balance model using
  greenhouse forcing, ice-albedo feedback, graph diffusion, and lapse-rate
  cooling.
- Four-season temperature now follows seasonal insolation, obliquity,
  eccentricity, land-ocean thermal inertia, and maritime lag.
- Winds now have a seasonal background circulation with migrating ITCZ and
  winter-strengthened storm tracks, plus bounded land-sea thermal, orographic,
  and C4c pressure/SST weak-coupled anomalies.  C5c2 adds a
  coupled-consistency gate so pressure, wind, moisture access, monsoon
  potential, and precipitation must support each other rather than pass as
  isolated fields.
- Ocean currents now use a reduced basin-streamfunction proxy that never
  crosses the solved ocean mask, produces gyres, boundary-current classes,
  strait exchange, SST anomalies, and current heat transport.  C5b1 adds an
  OSCAR/OISST spatial gate so earthlike strongest currents are checked for
  near-coast/boundary-current placement rather than only speed magnitude.
  C5e9 strengthens the bounded gyre/boundary-current heat anomaly and adds
  same-latitude SST residual checks so the reduced SST map cannot pass as a
  mostly latitude-banded texture.
- Precipitation is now solved as four seasonal fields and aggregated back to
  annual compatibility fields.  It uses routed moisture access, ITCZ/storm-track
  proxies, monsoon potential, and climate-scale orographic enhancement/rain
  shadow.  C5d1 adds a seasonal-hydro placement gate so wet/dry cells and
  wet-season timing must be explainable by these generated support fields.
  C5e1 adds a conservative precipitation-pressure/wind feedback pass before the
  final hydroclimate solve, so precipitation now participates in the reduced
  coupling loop instead of only consuming pressure/wind/moisture outputs.  C5e3
  adds a formal coupling-convergence gate and suppresses planet-scale hydro
  feedback from tiny island/waterworld land fractions.  C5e4 expands the
  one-shot feedback into a 3-pass bounded hydroclimate feedback loop with an
  archived iteration-delta diagnostic.  C5e5 adds a bounded evaporation-SST
  heat-flux feedback inside the ocean-atmosphere coupling loop.  C5e6 makes the
  wind-stress/current response explicit, archived, and gated.  C5e7 adds an
  explicit source-ocean-basin to receiver-catchment accounting layer so wet
  response regions can be checked against diagnosed seasonal basin supply
  rather than only local landward pathway strength.  C5e8 feeds that accounting
  back into precipitation as a second, bounded, land-only redistribution that
  preserves each seasonal local moisture-budget mean.
- Rendered assets include annual and seasonal temperature, seasonal SST,
  seasonal wind, pressure/anomaly diagnostics, currents, ocean heat transport,
  ocean heat flux, evaporation, upwelling, and annual precipitation.

Main observed problems:

- A shared geography primitive layer now exists for continents, coastlines,
  basin ids, shelves, straits, barriers, and wind gaps.  C4a pressure/moisture,
  C4b currents, C4c weak coupling, and the first C4d regional precipitation
  response now consume it.  C4d now also emits a first object layer for
  seasonal monsoon, storm-track, rain-shadow, and wet/dry response regions;
  the object-level, placement-proxy, map-readability, and rendered contact-sheet
  gates now pass.  The contact sheets were visually reviewed and are sufficient
  for the current Earth-fitting pass.  C4e adds explicit seasonal
  moisture-flow-network objects on top of C4d, and C4f now uses those routed
  networks as a conservative first-pass seasonal land-precipitation
  redistribution response.  C4g adds local moisture-budget region ids and
  preserves each season's precipitation mean inside those land budget regions,
  not only across global land.  C4h then safely splits only large, coherent
  moisture-flow networks into local halo sectors, while small or unstable
  networks fall back to the continent budget.  C4i adds dominant source-ocean
  basin attribution for seasonal moisture pathways and uses it as a guard on
  sector splitting.  C4j adds final wet/dry precipitation-response region
  objects that bind the active response to source basins, budget regions, and
  moisture-flow networks.  Ocean precipitation remains unchanged, and the
  accepted Earth-fitting gate suite stays green.  The dedicated moisture
  response gate now checks C4f response archive presence plus C4g/C4h/C4i
  local-budget/source-basin coherence and C4j precipitation-response object
  continuity, conservation, sector splitting, pathway/source/support coupling,
  map readability, and waterworld false-positive limits.
- C4b/C4c/C5b1 currents and SST are still reduced: they solve a basin
  streamfunction and bounded SST/wind feedback, and C5b1 now checks spatial
  boundary-current placement against broad OSCAR/OISST envelopes, but this is
  not a dynamic ocean model with salinity, thermohaline circulation, or monthly
  current observations.
- Monsoon potential is now geography-derived and consumes bounded SST/current
  feedback.  C4d emits first-pass monsoon rainfall corridors, storm-track
  rainfall corridors, rain-shadow index, and conservative regional precipitation
  response plus `climate.hydroclimate_regions` objects.  C4e emits seasonal
  moisture-source/pathway fields and `climate.moisture_flow_networks` objects,
  while C4f actively shapes seasonal precipitation through
  `climate.moisture_flow_precipitation_response`, C4g constrains that shaping
  with `climate.moisture_budget_region_id`, C4h splits stable large
  moisture-flow networks into budget sectors, C4i adds
  `atmosphere.moisture_source_basin_id` so those sectors carry a dominant
  source-ocean basin attribution, and C4j exposes final response patches as
  `climate.precipitation_response_regions`.  C5e7 now adds
  `climate.source_basin_supply_index` and
  `climate.receiver_catchment_supply_balance` so receiver catchments carry a
  bounded source-supply ledger.  This is still diagnostic, not a strict water
  conservation solve.  C5c2 tightens the monsoon moisture gate so high monsoon
  potential is not generated in dry interiors solely from thermal low pressure.
- Precipitation is solved as four seasonal fields, but still has residual
  terrain-line or broad-band artifacts because ITCZ, storm tracks, monsoon
  inflow, and orographic rain are not yet solved as a single regional field.
- Some extreme presets can end with a stale climate ocean mask relative to final
  terrain if post-climate terrain/ocean-mask drift does not cross the current
  climate re-solve trigger.
- Real-Earth R4 comparison originally showed a dominant hydroclimate blocker:
  earthlike terminal worlds had land precipitation around 20% of the Earth
  baseline, near-zero forest/tropical biome area, and excess desert area.  The
  F1-F5 fitting pass has cleared the scalar, pattern, coarse biome,
  spatial-biome, seasonal-subtype, mountain-zonation, and windward/leeward
  gates.  The C5e9 R6 acceptance bundle now also clears ocean-spatial,
  coupled-consistency, seasonal-hydro placement, hydro-region, moisture-flow,
  moisture-response, receiver-catchment, and coupling-convergence gates with
  zero failures, warnings, or skipped checks.  C5e7/C5e8 visual comparison previews
  are now backfilled from archived `terminal_climate_arrays.npz` when a replay
  was run without rendering, so the Earth-vs-generated contact sheet contains
  real Earth plus all six generated worlds.  C5e9 keeps that evidence path and
  makes the annual temperature/SST comparison less purely zonal without broad
  precipitation lifting.  Remaining work should focus on geography-coupled
  climate mechanics rather than scalar rain or biome retuning.

## Earth-Based Fitting Track

Detailed plan:

- `docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md`

This track starts from the R4/R6 Earth reference library.  Tectonics and terrain
remain frozen.  During R0-R8, the fitting target is not a generated world: it is
one real-Earth subgraph at a time, replayed on the real-Earth grid and compared
directly against the matching Earth reference map.  The six accepted terminal
worlds are held back for R9 generalization guardrails only.  The old F0-F5
fitting notes remain historical regression evidence, but new repairs must
follow the Replay-R0-R9 geoscience subgraph order archived in
`docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md`.

Authoritative replay order:

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

Ordering rule:

- Downstream products may be rendered and scored as observers, but they cannot
  be tuned until their upstream support fields pass.
- Every active phase must be accepted from map-level evidence first: choose one
  real-Earth subgraph, compare the real-Earth reference map with the
  corresponding real-Earth replay map on the same Earth grid, write down the
  visible spatial residuals, then assign those residuals to the active or
  upstream mechanism before changing code.  Global averages, percentile
  envelopes, correlations, and pass/fail tables are only regression aids.
- Generated-world maps are not used to choose R0-R8 fixes.  They are rerun only
  as R9 promotion guardrails after the relevant real-Earth subgraph is
  plausible.
- Do not run or tune Aevum virtual-world climate maps while an R0-R8 subgraph
  is being fitted.  The only maps that can justify a mechanism change in those
  phases are the real-Earth reference subgraph, the real-Earth replay subgraph,
  and their residual/error maps.
- Sea ice, snow, clouds, vegetation feedback, Koppen, and biomes are blocked
  until R2a pressure/source geometry, R2b wind translation, R3 currents, R4
  SST/energy, R5 moisture transport, and R6 precipitation are physically
  plausible on real-Earth replay.

Current active fitting packet:

- Active subgraph: R2a seasonal SLP / pressure-source geometry replay on real
  Earth.
- Reference: `earth__seasonal_slp_anomaly_hPa` from the R6 Earth reference
  package.
- Replay: `atmosphere__land_sea_pressure_proxy` /
  `atmosphere__seasonal_pressure_proxy` from `real-earth-climate-replay`.
- Required visual evidence: Earth seasonal standardized SLP anomaly, replay
  seasonal standardized pressure-proxy anomaly, standardized pressure residual,
  pressure zonal-anomaly maps, and pressure-center / stationary-wave support
  maps when available.
- Current known residual: the R2a baseline is
  `out_real_earth_pressure_replay_r2a_current_20260706/`.  Replay pressure
  captures broad seasonal land thermal contrast, but it is too smooth and too
  blocky.  It lacks enough ocean basin pressure centers and mountain/coast
  stationary-wave structure to organize later storm tracks.  The previous best
  wind checkpoint,
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/`, remains a
  downstream observer and historical reason for moving upstream into pressure;
  it is not the current fitting target.
- Current code owner: R2a pressure-source geometry.  The next mechanism should
  first expose explicit basin-scale pressure-center / stationary-wave objects,
  then judge any pressure-source change by the pressure contact sheet.  R2b
  wind translation is blocked until R2a pressure maps stop reading as
  over-smooth continent/ocean blobs.
- Blocked next phases: R2b wind translation, R3 OSCAR currents, R4 OISST SST,
  R5 moisture, R6 precipitation, R7 cryosphere/cloud/vegetation, and R8
  classes/biomes.

Historical first-pass R4 baseline conclusion:

- Temperature is not the largest blocker, although ocean temperatures should be
  monitored during later tuning.
- Current vectors are somewhat fast relative to the Earth drifter baseline, but
  their thermal effect is weak.
- The severe mismatch is land precipitation: `earthlike_seed42` and
  `earthlike_seed909` produce about 143-158 mm/yr mean land precipitation versus
  the Earth baseline of about 730 mm/yr.
- Forest and tropical biome tuning is blocked until the water-cycle fields are
  plausible.

F0 implementation status:

- `aevum.diagnostics.earth_climate_fitting` writes the phase-priority report,
  run CSV, candidate lever CSV, guardrail CSV, and cross-preset guardrail
  verdict.
- CLI command: `aevum earth-climate-fit-report`.
  Use `--fail-on-guardrail` for automated sweeps; warnings remain non-fatal.
- `aevum.diagnostics.terminal_climate_replay` reruns climate and static biomes
  on frozen terminal arrays without rerunning tectonics or terrain.
- CLI command: `aevum terminal-climate-replay`.

F1-F5 first fitting pass, 2026-07-05:

- Climate and biome edits were made only in `aevum.modules.climate` and
  `aevum.modules.biosphere`; accepted plate/terrain generation was not changed.
- F1: ocean cells now use a mixed-layer/SST floor near seawater freezing for
  calibration against OISST, while sea-ice state still comes from the raw cold
  climate solution.  Earthlike global temperature deltas improve to about
  -0.21 C and -1.69 C for the two current seeds.
- F2: exported near-surface ocean-current vectors are calibrated separately from
  the reduced heat-transport proxy.  C5b1 uses the R6 OSCAR/OISST baseline and
  now checks both p90 speed magnitude and near-coast/far-ocean strong-current
  placement; earthlike current p90 ratio is about `0.90`.
- F3: moisture access now includes a weak free-atmosphere reservoir and stronger
  source-ocean/monsoon coupling.  Earthlike land moisture-access p75 is now
  about 0.59-0.65.
- F4: land precipitation now uses stronger routed moisture rainout, nonlinear
  ITCZ/storm-track/monsoon wet-core enhancement, a warm high-access convective
  tail, and stronger geography-conditioned subtropical subsidence.  Earthlike
  dry flags are cleared; Pattern7 global land precipitation mean/p50/p90
  ratios are about `0.57/0.62/0.57` and `0.69/0.84/0.70`.  C5d1 now verifies
  seasonal hydroclimate placement against the generated moisture, ITCZ,
  monsoon, storm-track, rain-shadow, and response fields.
- F5: biome thresholds are less conservative, generalization preserves
  climate-supported forest/tropical semantic patches, cold-dry land is no
  longer overwritten as desert, and moist-temperate forest now begins at
  520 mm/yr.  The coarse Earth biome envelope gate passes against both Koppen
  proxy and RESOLVE references.  A follow-up spatial-biome gate checks latitude
  organization; cool climates now use lower forest/desert precipitation
  thresholds and high-latitude cold-dry land is classified as tundra/ice before
  desert.  A seasonal-subtype gate now checks Koppen-like dry/wet-quarter
  organization; low-latitude seasonal contrast is strengthened without changing
  annual precipitation.  A mountain-zonation gate now checks high-mountain
  alpine ecology and desert excess; cool high-elevation land enters alpine/
  tundra semantics before arid classification.  A windward/leeward gate now
  checks seasonal wind against mountain-slope precipitation; climate now uses
  slope-wind exposure to redistribute seasonal land precipitation from leeward
  to windward slopes while preserving the seasonal land mean.
- Current first-pass report:
  `out_earth_climate_fitting_f1_oceanfloor_20260705/earth_climate_fitting_report.md`.
- Current guardrail report:
  `out_earth_climate_fitting_f5wind3_gate_20260705/earth_climate_fitting_report.md`.
- Current guardrail CSV:
  `out_earth_climate_fitting_f5wind3_gate_20260705/earth_climate_guardrails.csv`.
- Current Earth-pattern gate report:
  `out_earth_climate_pattern_gate_f5wind3_20260705/earth_climate_pattern_gate_report.md`.
- Current Earth-biome gate report:
  `out_earth_climate_biome_gate_f5wind3_20260705/earth_climate_biome_gate_report.md`.
- Current Earth-spatial-biome gate report:
  `out_earth_climate_spatial_biome_gate_f5wind3_20260705/earth_climate_spatial_biome_gate_report.md`.
- Current Earth-seasonal-subtype gate report:
  `out_earth_climate_seasonal_subtype_gate_f5wind3_20260705/earth_climate_seasonal_subtype_gate_report.md`.
- Current Earth-mountain-zonation gate report:
  `out_earth_climate_mountain_zonation_gate_f5wind3_20260705/earth_climate_mountain_zonation_gate_report.md`.
- Current Earth-windward/leeward gate report:
  `out_earth_climate_windward_leeward_gate_f5wind3_20260705/earth_climate_windward_leeward_gate_report.md`.
- Current rendered six-world climate assets:
  `out_terminal_climate_replay_f5wind3_render_20260705/`.

Residual fitting priorities:

- F1 temperature remains only a watch item because land/ocean meaning differs
  between generated worlds and the mixed WorldClim/OISST Earth reference.
- F4 climate-pattern failures are cleared at the current pattern-gate
  strictness: wet tropics, dry subtropics, high-latitude cold envelope, and
  mountain wet-tail checks now pass.
- F4 seasonal-hydro placement now passes on the C5c2 replay under the C5d1
  gate: wet cells have process support, dry cells have moisture/rain-shadow
  explanations, and wet seasons align with maximum support.  This is still a
  bounded weak-coupling solution rather than a full atmosphere-ocean GCM.
- F5 coarse biome-envelope checks now pass: earthlike forest+tropical land
  fractions are about 0.281 and 0.350 after spatial-biome tuning, clearing the
  gate floor against Koppen proxy and RESOLVE coarse references.
- F5 spatial-biome checks now pass: tropical biome is low-latitude, subtropical
  dry belts remain present, cool-midlatitude forest/desert fractions are no
  longer inverted, and high-latitude desert is near zero in both earthlike
  seeds.
- F5 seasonal-subtype checks now pass: low-tropical dry-quarter subtype area is
  present, low/mid-latitude precipitation seasonality stays within Earth
  tolerance, and tropical biome semantics remain in bounds.
- F5 mountain-zonation checks now pass: high-mountain alpine ecology is present
  in both earthlike seeds, high-mountain desert fractions fall below the
  failure threshold, and mountain cooling remains in bounds.
- Windward/leeward checks now pass: earthlike mountain-slope windward sides are
  wetter than leeward sides in annual and seasonal metrics without changing the
  accepted terminal terrain.
- Cross-preset scalar guardrails, Earth-pattern gate, Earth-biome gate,
  Earth-spatial-biome gate, Earth-seasonal-subtype gate, and
  Earth-mountain-zonation gate, Earth-windward/leeward gate, and
  Earth-seasonal-hydro placement gate currently pass with zero warnings and
  zero failures.

Current Earth-fitting acceptance:

- The current authoritative acceptance bundle is C5e9/R6 and uses the
  ocean-structure replay:
  `out_terminal_climate_replay_c5e9_ocean_structure_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e9_ocean_structure_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e9_ocean_structure_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.  All F1-F5
  phases are low-priority `watch`.
- The complete C5e9 acceptance suite passes with `0` failures, `0` warnings,
  and `0` skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement, hydro-region,
  coupling-convergence, moisture-flow, moisture-response, and
  receiver-catchment.
- C5e9 supersedes C5e8 by preserving the receiver-supply precipitation feedback
  while adding anti-zonal SST/heat-transport structure checks.  The next active
  climate-system target should build on this baseline rather than broadly
  lifting precipitation.

## Design Goals

The climate layer should remain reduced and fast, but it must represent the
major causal controls visible in planetary-scale maps:

1. Seasonal insolation from latitude, obliquity, eccentricity, and rotation mode.
2. Land-ocean heat capacity contrast and continentality.
3. Atmospheric circulation with seasonally migrating ITCZ, subtropical dry belts,
   mid-latitude storm tracks, and tidally locked variants.
4. Ocean heat transport by basin-aware current proxies.
5. Monsoon and coastal precipitation driven by seasonal pressure/thermal contrast.
6. Orographic rain and rain-shadow as regional effects, not single-cell stripes.
7. Ice/snow/cloud/vegetation feedbacks strong enough to soften hard thermal walls.
8. Assets that expose both annual means and seasonal structure.

The next architecture target is geography-coupled climate: ocean currents,
monsoon circulation, moisture transport, and precipitation must all consume the
same explicit geography diagnostics.  The intended dependency direction is:

```text
terrain + ocean mask
  -> continents / ocean basins / shelves / straits / coasts / barriers / passes
  -> seasonal pressure + moisture access + wind corridors
  -> basin gyres + SST/upwelling + monsoon inflow
  -> seasonal precipitation + runoff + biomes
```

Non-goals for the near term:

- Full 3D GCM dynamics.
- Cloud microphysics.
- ENSO/PDO/NAO-style internal variability.
- Detailed thermohaline circulation with salinity conservation.
- Daily weather, cyclones, or synoptic storm tracks.

## Target Climate Fields

Keep existing annual fields for downstream compatibility:

- `climate.surface_temperature`: annual mean K.
- `climate.precipitation`: annual total or annual-mean-equivalent mm/yr.
- `climate.evaporation`
- `climate.runoff`
- `atmosphere.wind`
- `ocean.currents`
- `cryosphere.sea_ice`
- `cryosphere.ice_sheet`

Add new seasonal/diagnostic fields:

- `climate.seasonal_temperature`: shape `(4, n_cells)`, K.
- `climate.seasonal_precipitation`: shape `(4, n_cells)`, mm/yr equivalent.
- `climate.temperature_seasonality`: annual max-min K.
- `climate.precipitation_seasonality`: max seasonal precip / annual mean.
- `atmosphere.seasonal_wind`: shape `(4, n_cells, 3)`.
- `atmosphere.background_seasonal_wind`: astronomical/latitude circulation,
  shape `(4, n_cells, 3)`.
- `atmosphere.thermal_wind_anomaly`: land-ocean thermal pressure anomaly,
  shape `(4, n_cells, 3)`.
- `atmosphere.orographic_wind_anomaly`: terrain steering/blocking anomaly,
  shape `(4, n_cells, 3)`.
- `atmosphere.land_sea_pressure_proxy`: seasonal thermal pressure proxy,
  shape `(4, n_cells)`.
- `atmosphere.geographic_circulation_index`: combined geography-driven
  circulation strength, shape `(n_cells,)`.
- `ocean.current_heat_transport`: scalar or vector diagnostic.
- `climate.continentality`: land-ocean thermal inertia / distance proxy.
- `climate.monsoon_index`: signed seasonal wet/dry contrast.

Planned shared geography and coupled-climate fields:

- `climate.continent_id`: connected exposed-land component id.
- `climate.continent_interiority`: distance/insulation from ocean, scaled by
  component size.
- `climate.coast_orientation`: tangent direction from land toward adjacent
  ocean and east/west/north/south-facing coast diagnostics.
- `ocean.basin_id`: connected ocean basin id, already present as a C3 field.
- `ocean.shelf_index`: broad shelf/slope/deep-ocean diagnostic from bathymetry
  and graph distance from land.
- `ocean.strait_index`: narrow ocean-gateway and exchange bottleneck diagnostic.
- `ocean.gyre_id`: basin gyre cell id.
- `ocean.current_streamfunction`: reduced basin streamfunction for current
  direction and speed.
- `ocean.sst_anomaly`: ocean-current/upwelling sea-surface-temperature anomaly.
- `terrain.barrier_index`: mountain/highland barrier strength for winds and
  moisture.
- `terrain.wind_gap_index`: low passes and corridors through barriers.
- `atmosphere.seasonal_pressure_proxy`: final pressure field after land-sea,
  SST, and terrain effects.
- `atmosphere.moisture_access`: seasonal ability of moist air to reach each
  cell from source oceans.
- `atmosphere.monsoon_potential`: seasonal geography-derived monsoon potential.
- `climate.seasonal_precipitation`: seasonal rainfall/snowfall, replacing the
  annual-only precipitation calculation as the primary hydroclimate solve.

Suggested season order:

1. DJF
2. MAM
3. JJA
4. SON

## Geography-Coupled Climate Architecture

The post-C3 work should not add independent rules for monsoons, currents, and
rainfall.  It should introduce a shared geography diagnostic layer that all
later climate components consume.  This keeps causality inspectable and reduces
the risk that separate heuristics disagree about the same coastline, basin, or
mountain belt.

Core design constraints:

- All geography diagnostics operate on the sphere graph and must be continuous
  across the dateline.
- Basin, continent, coastline, shelf, strait, barrier, and pass objects are
  computed once per relevant terrain/ocean state and reused by atmosphere,
  ocean, hydroclimate, biosphere, and map compiler.
- Seasonal pressure, wind, currents, SST anomaly, evaporation, and precipitation
  are solved as a weakly coupled reduced system.  They do not need full GCM
  dynamics, but they must exchange the right first-order information.
- Energy and water budgets remain bounded.  Ocean heat transport redistributes
  heat; it should not add net global energy.  Seasonal precipitation should
  aggregate back to the annual fields expected by downstream modules.
- Every derived map should be explainable: a cell should be dry because it is
  far from source oceans, behind a barrier, under a cold current, or outside the
  seasonal convergence path, not because of an opaque local color rule.

Shared diagnostic objects:

- Continents: connected land components with area, centroid, latitudinal span,
  interiority, coastline length, adjacent ocean basins, and seasonal thermal
  inertia.
- Ocean basins: connected ocean components with area, open/closed status,
  constrictions, shelves, deep interiors, polar/tropical connections, and
  adjacent continents.
- Coasts: paired land-ocean boundary cells with tangent orientation, oceanward
  normal, alongshore direction, facing classification, shelf width, and adjacent
  basin id.
- Straits and gateways: narrow ocean corridors where current exchange and
  moisture access should be limited or focused.
- Terrain barriers: smoothed mountain/highland belts with barrier strength,
  cross-barrier normal, along-barrier direction, and low-pass/gap diagnostics.
- Moisture corridors: graph paths from warm source oceans into land, damped by
  barriers and boosted by seasonal onshore flow.

Weak coupling loop:

```text
seasonal insolation + geography primitives
  -> land/ocean/SST pressure proxy
  -> terrain-steered seasonal wind
  -> basin streamfunction + upwelling + SST anomaly
  -> evaporation + source-ocean humidity
  -> moisture access + convergence + seasonal precipitation
  -> annual temperature/precip/runoff compatibility fields
```

Implementation should use 2-4 bounded iterations for this loop after all
diagnostic fields are present.  Each iteration should preserve tangent vectors,
avoid crossing land barriers, and cap anomalies before they feed the next field.

## Implementation Plan

### Phase C0 - Baseline Diagnostics

Status: complete

Purpose: make current climate shortcomings measurable before changing the model.

Tasks:

- Add diagnostics for latitudinal temperature jumps, land/ocean seasonal
  amplitude placeholder, and precipitation striping.
- Add helper metrics in `validation.py`:
  - `max_adjacent_lat_band_delta_C`
  - `land_ocean_temperature_contrast_C`
  - `precip_orographic_concentration`
  - `coastal_temperature_asymmetry_index`
- Add render support for optional diagnostic maps if fields exist.

Tests:

- Earthlike annual mean temperature should not have adjacent 10-degree band
  jumps above a configured threshold.
- Precipitation should not concentrate too much of land rainfall in narrow
  orographic bands.

Assets to inspect:

- `temperature.png`
- `precip.png`
- `biomes.png`

Risks:

- Diagnostics may initially fail on current output.  Treat them as warnings
  until phases C1-C4 land.

Implementation notes:

- `validation.climate_diagnostics()` now reports temperature structure,
  precipitation/orographic concentration, coastal temperature asymmetry, and
  seasonal-field availability.
- `validation.check_climate_diagnostics()` treats invalid data as hard failures
  and current model-quality gaps as warnings.
- `render_world()` can render optional scalar climate diagnostics when future
  phases add them to the world fields.

### Phase C1 - Seasonal Insolation and Thermal Inertia

Status: complete

Purpose: replace the single annual thermal state with a four-season climatology.

Tasks:

- Implement seasonal solar factor from latitude and obliquity.
- Include eccentricity as a small seasonal amplitude modifier.
- Add land/ocean effective heat capacity:
  - ocean: high inertia, low seasonal amplitude.
  - land: lower inertia, larger seasonal amplitude.
  - ice: high albedo and seasonal persistence.
- Solve four seasonal temperature snapshots, then aggregate annual mean.
- Keep annual `climate.surface_temperature` as the average of the seasonal field.

Tests:

- For Earthlike, northern land JJA mean should exceed northern land DJF mean.
- Land seasonal amplitude should be greater than ocean seasonal amplitude.
- Annual mean temperature remains within existing plausibility bounds.
- Tidally locked preset should still use day/night forcing rather than normal
  Earthlike seasons.

Assets:

- Add `temperature_seasons.png`.
- Add `temperature_seasonality.png`.

Expected visual result:

- The climate map should stop reading as a single static hot belt.
- High latitudes should show seasonal contrast rather than a hard permanent wall.

Implementation notes:

- `climate.seasonal_temperature` now stores DJF/MAM/JJA/SON temperature fields
  with shape `(4, n_cells)`.
- `climate.surface_temperature` is now the mean of the four seasonal fields.
- `climate.temperature_seasonality` stores per-cell max-min seasonal amplitude.
- `climate.continentality` stores the land-ocean thermal inertia proxy.
- Tidally locked worlds keep a static day/night field rather than Earthlike
  obliquity seasons.
- `render_world()` now writes `temperature_seasons.png` when seasonal
  temperature is available.

### Phase C2 - Seasonal Winds, ITCZ, and Storm Tracks

Status: complete

Purpose: make atmospheric circulation move with the seasons.

Tasks:

- Replace fixed wind bands with season-dependent winds:
  - ITCZ follows seasonal solar maximum.
  - Hadley cell edge shifts by season and rotation rate.
  - Subtropical dry belts shift seasonally.
  - Mid-latitude storm tracks strengthen in winter hemisphere.
- Keep simple vector winds in tangent coordinates.
- Add `atmosphere.seasonal_wind`; annual `atmosphere.wind` remains a mean.

Tests:

- ITCZ/convergence intensity shifts north in JJA and south in DJF for non-locked
  worlds with Earthlike obliquity.
- Mid-latitude winter storm-track intensity is stronger than summer counterpart.
- Wind vectors remain tangent and finite.

Assets:

- Add `wind_seasons.png` or sparse quiver overlays.
- Add optional `itcz_track.png`.

Expected visual result:

- Tropical rain belt becomes seasonally mobile.
- Subtropical dry zones become less fixed and less wall-like.

Implementation notes:

- `atmosphere.seasonal_wind` now stores DJF/MAM/JJA/SON wind vectors with shape
  `(4, n_cells, 3)`.
- `atmosphere.wind` is now the four-season mean for downstream compatibility.
- `atmosphere.itcz_latitude` stores the four seasonal ITCZ latitudes.
- `atmosphere.itcz_intensity` and `atmosphere.storm_track_intensity` provide
  seasonal convergence/storm-track drivers for C4 precipitation.
- Tidally locked worlds keep a static day/night circulation rather than
  Earthlike seasonal wind migration.

### Phase C2.5 - Geography-Driven Circulation Anomalies

Status: complete

Purpose: make the C2 background circulation respond to the actual geography.
The current C2 wind/ITCZ/storm-track fields are astronomical and zonal; this
phase adds land-ocean thermal pressure, terrain steering, and basin/coast
structure before C3/C4 consume the circulation fields.

Design principle:

- Keep C2 as the large-scale planetary background.
- Add bounded regional anomalies from geography.
- Preserve tangent vector fields and downstream compatibility:
  `atmosphere.wind` remains annual mean and `atmosphere.seasonal_wind` remains
  the final four-season wind.

Inputs:

- `climate.seasonal_temperature`
- `climate.continentality`
- `terrain.elevation_m`
- `ocean.mask`
- `terrain.basins` if available
- `atmosphere.background_seasonal_wind` from C2
- Later C3 input: `ocean.current_heat_transport`

Outputs:

- `atmosphere.background_seasonal_wind`: copy of the C2 astronomical wind.
- `atmosphere.land_sea_pressure_proxy`: seasonal thermal low/high proxy.
- `atmosphere.thermal_wind_anomaly`: inflow/outflow caused by land-ocean
  thermal contrast.
- `atmosphere.orographic_wind_anomaly`: terrain steering and barrier damping.
- `atmosphere.geographic_circulation_index`: annual diagnostic strength of
  geography-driven circulation.
- Updated `atmosphere.seasonal_wind`: background plus bounded anomalies.

Algorithm:

1. Compute seasonal thermal pressure proxy:
   - Start from seasonal temperature anomaly relative to each cell's annual mean.
   - Amplify land by `continentality`; damp open ocean.
   - Smooth on the graph to create regional thermal lows/highs.
   - Summer continental heat lows should be negative pressure proxy; winter cold
     continents should be positive pressure proxy.
2. Derive thermal wind anomaly:
   - Compute graph-gradient flow from high pressure to low pressure.
   - Project gradients onto local tangent vectors.
   - Limit anomaly magnitude to a fraction of background wind speed.
   - Apply only where land/ocean contrast is meaningful; waterworld should stay
     weak.
3. Detect coastal inflow/outflow:
   - Identify land cells adjacent to ocean and ocean cells adjacent to land.
   - Summer: add ocean-to-land inflow toward continental thermal lows.
   - Winter: add land-to-ocean outflow from continental cold highs.
   - Scale by continent size/continentality and nearby warm-ocean availability.
4. Add terrain steering:
   - Use smoothed elevation gradients, not raw cell relief.
   - Dampen wind component crossing major mountain belts.
   - Add along-barrier deflection so air routes around ranges rather than
     punching through as a stripe.
5. Compose final wind:
   - `seasonal_wind = background + thermal_anomaly + orographic_anomaly`.
   - Re-project to tangent plane.
   - Smooth lightly and cap magnitude.
   - Annual `atmosphere.wind` remains seasonal mean.

Tests:

- Seasonal wind vectors are tangent and finite.
- `atmosphere.seasonal_wind.mean(axis=0)` equals `atmosphere.wind`.
- Large summer-heated continents have stronger onshore flow than winter at
  adjacent coasts.
- Winter continents produce net offshore anomaly relative to summer.
- Waterworld has weak `geographic_circulation_index` and weak thermal anomaly.
- Arid/large-continent worlds show stronger geographic circulation than
  waterworld.
- Orographic anomaly reduces cross-mountain wind component without making wind
  zero everywhere.

Assets:

- Add `thermal_wind_anomaly.png` or seasonal panel.
- Add `land_sea_pressure.png`.
- Add `geographic_circulation_index.png`.

Expected visual result:

- Same latitude regions no longer share identical wind direction/intensity.
- Continents create summer inflow and winter outflow.
- Mountain belts steer winds regionally instead of only creating precipitation
  stripes later.
- Waterworld retains mostly zonal background circulation.

Implementation notes:

- C2 background winds are now preserved in
  `atmosphere.background_seasonal_wind`.
- `atmosphere.land_sea_pressure_proxy` is derived from seasonal temperature
  anomaly and continentality, with fixed-scale normalization so waterworlds do
  not get artificial monsoon-strength anomalies.
- `atmosphere.thermal_wind_anomaly` adds bounded high-to-low pressure flow plus
  smoothed coastal inflow/outflow.
- `atmosphere.orographic_wind_anomaly` damps cross-barrier wind and adds limited
  along-barrier steering from smoothed topography.
- `atmosphere.seasonal_wind` is the final composed wind; `atmosphere.wind`
  remains its annual mean.

### Phase C3 - Ocean Current Heat Transport

Status: complete

Purpose: make oceans actively shape temperature, especially coasts.

Dependency note: C3 should consume the final C2.5 seasonal wind where possible,
not the purely zonal C2 background.  This lets geography-driven monsoon/coastal
winds shape gyres, upwelling, and current heat transport.

Tasks:

- Build ocean basin connectivity from `ocean.mask`.
- Generate wind-driven gyre proxies:
  - low-latitude westward flow.
  - mid-latitude eastward return flow.
  - western boundary intensification heuristic.
  - coastal upwelling on appropriate eastern boundaries.
- Add polarward heat transport along western boundary currents.
- Add equatorward cooling along eastern boundary currents.
- Feed current heat transport into seasonal temperature iteration.

Tests:

- Currents never cross land.
- Ocean heat transport is finite and mostly follows ocean-connected cells.
- Earthlike-style worlds show nonzero east/west coast thermal asymmetry.
- Coastal cold-current zones should reduce evaporation and rainfall locally.

Assets:

- Add `currents.png`.
- Add `ocean_heat_transport.png`.

Expected visual result:

- Same-latitude coasts no longer have identical temperatures.
- Some western ocean margins become cooler/drier; some eastern continental
  margins become warmer/wetter depending on hemisphere and flow.

Implementation notes:

- `ocean.basin_id` now records connected ocean-basin components on the sphere
  graph; `ocean.solved_mask` records the ocean mask used by the latest climate
  solve so diagnostics can separate C3 topology errors from later terrain drift.
- `ocean.currents` is now a basin-constrained tangent vector field, not a scaled
  copy of annual wind.  It combines low-latitude westward drift, mid-latitude
  eastward return flow, boundary-current intensification, and bounded smoothing
  only within connected ocean cells.
- `ocean.current_heat_transport` stores the current-induced temperature anomaly
  proxy, with ocean heat redistributed into adjacent coasts and global mean
  removed to avoid adding net energy.
- `ocean.upwelling` marks eastern-boundary cold-current/upwelling zones.
- Seasonal temperature consumes `ocean.current_heat_transport`; annual
  `climate.surface_temperature` remains the four-season mean.
- Hydroclimate now applies a modest cold-current drying adjustment to
  evaporation and precipitation, while full seasonal hydroclimate remains C4.
- `render_world()` now writes `currents.png`, `ocean_heat_transport.png`,
  `upwelling.png`, and `ocean_basin_id.png` when those fields exist.
- `validation.climate_diagnostics()` now reports ocean-current shape, tangent,
  solved-mask crossing, final terrain mismatch, heat-transport, upwelling, and
  basin metrics.

### Phase C3.5 - Shared Geography Primitive Layer

Status: complete

Purpose: build the common geography substrate that ocean currents, seasonal
pressure, monsoons, moisture transport, precipitation, biomes, and the compiler
all consume.  This phase is intentionally before new hydroclimate work; without
it, later improvements will keep adding separate local heuristics.

Inputs:

- `terrain.elevation_m`
- `ocean.sea_level_m`
- `ocean.mask`
- `crust.type`
- `terrain.basins` if available
- Existing graph topology and cell areas

Outputs:

- `climate.continent_id`
- `climate.continent_interiority`
- `climate.coast_orientation`
- `climate.coast_distance`
- `ocean.basin_id` refinement
- `ocean.shelf_index`
- `ocean.strait_index`
- `terrain.barrier_index`
- `terrain.wind_gap_index`
- Object sets for continents, ocean basins, coastline segments, straits, and
  major barrier belts.

Algorithm:

1. Compute connected land and ocean components on the sphere graph.
2. Create continent objects:
   - area fraction, centroid, latitudinal span, compactness, coastal perimeter,
     adjacent ocean basins, and interior cells.
   - interiority from graph distance to ocean, scaled by continent area and
     reduced by lakes/seaways if present.
3. Create ocean basin objects:
   - area fraction, latitudinal span, open/closed status, shelf fraction, deep
     interior, adjacent continents, and connection graph through straits.
4. Classify coasts:
   - oceanward normal, alongshore direction, east/west/north/south facing,
     adjacent basin id, local shelf width, and whether coast is exposed to open
     basin or enclosed sea.
5. Detect shelves and straits:
   - shelf from depth and graph distance to land.
   - strait from narrow ocean corridors where deleting a small cell set greatly
     changes basin connectivity or local ocean width is low.
6. Detect terrain barriers and wind gaps:
   - smooth elevation to climate scale.
   - derive barrier strength from high, continuous relief belts.
   - identify low passes/gaps where moisture/wind can cross.

Tests:

- Land and ocean components cover all cells exactly once and are dateline
  continuous.
- `ocean.basin_id` never labels land as ocean.
- Coast orientation vectors are tangent and point from land toward adjacent
  ocean.
- Shelf index is highest near coasts and lower in deep basin interiors.
- Strait index is nonzero only in narrow ocean corridors and remains zero on
  land.
- Barrier index is concentrated on broad highlands/mountain belts, not one-cell
  spikes.
- Waterworld has near-zero continent objects but valid ocean-basin diagnostics.
- Arid/large-continent worlds produce strong interiority.

Assets:

- `geography_primitives.png`
- `continent_id.png`
- `ocean_basins.png`
- `shelf_strait.png`
- `coast_orientation.png`
- `terrain_barriers.png`
- `wind_gaps.png`

Expected visual result:

- The maps reveal interpretable continents, open basins, marginal seas,
  shelves, straits, and major barriers.
- Later climate maps can be audited against these primitives.

Implementation notes:

- `ClimateModule` now computes shared geography primitives during each climate
  solve.  These are diagnostic fields at C3.5; they do not yet change the
  climate solution.
- Added land/ocean connected components, normalized distance-to-coast,
  continent interiority, coastal orientation/strength, shelf index, strait
  index, terrain barrier index, and wind-gap index.
- Added object summaries for `climate.continents`, `ocean.basins`,
  `climate.coastline_segments`, `ocean.straits`, and
  `terrain.barrier_belts`.
- `render_world()` now writes `geography_primitives.png`, `continent_id.png`,
  `continent_interiority.png`, `coast_distance.png`, `coast_strength.png`,
  `coast_orientation.png`, `shelf_index.png`, `strait_index.png`,
  `terrain_barriers.png`, and `wind_gaps.png`.
- `validation.climate_diagnostics()` now reports a `geography` section covering
  component coverage, invalid shapes, nonfinite values, coast-vector tangency,
  shelf/deep-ocean contrast, strait land leakage, barrier/highland contrast, and
  object counts.
- Strait detection was iterated once after visual review: isolated one-cell
  coastal noise is filtered out, leaving fewer high-confidence narrow-gateway
  candidates.  The field remains a proxy and should be replaced/refined by C4b
  streamfunction/gateway logic.

### Phase C4a - Geography-Derived Seasonal Pressure and Moisture Access

Status: first pass complete; Earth/preset gate passing

Purpose: replace the current local land-sea pressure anomaly with a pressure and
moisture-access system derived from continent objects, nearby source oceans, SST
anomaly, terrain barriers, and seasonal heating.

Dependency note: C4a consumes C3.5 primitives and existing C1/C2/C2.5 seasonal
temperature/wind.  It should preserve C2 as the planetary background but make
regional pressure and moisture access explicitly geography-derived.

Outputs:

- `atmosphere.seasonal_pressure_proxy`
- Updated `atmosphere.land_sea_pressure_proxy` or migration path toward replacing
  it.
- `atmosphere.moisture_access`
- `atmosphere.monsoon_potential`
- `atmosphere.source_ocean_warmth`
- `atmosphere.terrain_blocking`

Algorithm:

1. Compute seasonal thermal pressure over land:
   - seasonal temperature anomaly, continent interiority, land heat capacity,
     and component size produce summer heat lows and winter cold highs.
2. Compute ocean pressure/source modifiers:
   - nearby warm SST/upwelling/cold-current anomalies change evaporation and
     low-level pressure over source oceans.
   - enclosed seas have limited source strength unless warm and connected to
     large open basins.
3. Route moisture access:
   - start from warm ocean/source cells.
   - propagate over graph along seasonal wind direction.
   - damp across `terrain.barrier_index`, restore through `wind_gap_index`, and
     decay with inland distance.
4. Diagnose monsoon potential:
   - positive where summer continental heat low, onshore pressure gradient,
     warm adjacent source ocean, and passable terrain corridor align.
   - negative/offshore where winter high pressure pushes flow away from land.
5. Compose updated regional wind anomaly:
   - gradient flow from high to low pressure.
   - terrain steering around barriers and through gaps.
   - cap and smooth anomalies before combining with background wind.

Tests:

- For Earthlike-like worlds, major summer continents show negative pressure
  proxy and higher monsoon potential near warm adjacent basins.
- Winter versions of the same continents show offshore or reduced monsoon
  potential.
- Waterworld has near-zero monsoon potential and weak continent-driven pressure.
- Arid/large-continent worlds have strong interior pressure contrast but limited
  moisture access where source oceans are far or barriers block flow.
- Moisture access decreases across high barriers and recovers through wind gaps.
- Pressure and wind anomaly vectors remain finite, tangent, and dateline
  continuous.

Implementation notes:

- `ClimateModule` now emits `atmosphere.seasonal_pressure_proxy`,
  `atmosphere.moisture_access`, `atmosphere.monsoon_potential`,
  `atmosphere.source_ocean_warmth`, and `atmosphere.terrain_blocking`.
- `aevum.diagnostics.earth_climate_monsoon_moisture_gate` checks the C4a
  fields against broad Earth monsoon envelopes and preset guardrails.
- The gate caught a real failure in `f5wind3`: `waterworld_seed707` produced
  continent-scale monsoon potential and pressure reversal on small islands.
- The fix adds absolute continent-area damping to continent interiority and
  seasonal pressure normalization, so tiny archipelagos cannot be promoted to
  continental heat lows merely because they are the largest land in a
  waterworld.
- Current replay:
  `out_terminal_climate_replay_c4a1_20260705/`.
- Current gate:
  `out_earth_climate_monsoon_moisture_gate_c4a1_20260705/`, verdict `pass`
  with `0` failures and `0` warnings.
- Key C4a metrics: earthlike summer monsoon-potential p90 is `0.359` and
  `0.679`; waterworld summer monsoon-potential p90 is now `0.059` and `0.113`;
  waterworld seasonal pressure reversal is now weak at about `-0.087` and
  `-0.217`.

Assets:

- `seasonal_pressure.png`
- `moisture_access.png`
- `monsoon_potential.png`
- `source_ocean_warmth.png`
- `terrain_blocking.png`

Expected visual result:

- Monsoon-capable coasts are visible before precipitation is solved.
- Large continents, warm adjacent seas, and mountain barriers explain where
  seasonal inflow can and cannot reach.

### Phase C4b - Basin-Streamfunction Ocean Currents

Status: first pass complete; six-world gate regression passing

Purpose: replace the C3 boundary-current heuristic with basin-scale current
cells derived from ocean-basin geometry, wind-stress curl, boundaries, shelves,
and strait exchange.

Dependency note: C4b should keep the C3 fields for compatibility but change how
they are generated.  C3 remains the fallback proxy while this phase is developed.

Outputs:

- `ocean.gyre_id`
- `ocean.current_streamfunction`
- Updated `ocean.currents`
- Updated `ocean.current_heat_transport`
- Updated `ocean.upwelling`
- `ocean.boundary_current_type`
- `ocean.strait_exchange`
- `ocean.sst_anomaly`

Algorithm:

1. For each open ocean basin, compute wind-stress curl from seasonal or annual
   wind projected onto local east/north tangent basis.
2. Solve a reduced graph Poisson/relaxation problem for streamfunction within
   each basin:
   - land boundaries are no-flow.
   - straits allow limited exchange weighted by `ocean.strait_index`.
   - shelves slow currents and focus alongshore flow.
3. Derive current vectors from streamfunction gradients and re-project tangent.
4. Add western-boundary intensification from basin geometry:
   - identify basin western boundaries from basin centroid and planetary
     rotation frame, not only local coast orientation.
   - intensify narrow boundary jets where streamfunction return flow hits a
     coast.
5. Diagnose eastern-boundary cold currents and upwelling:
   - require alongshore equatorward flow, appropriate wind direction/curl,
     subtropical latitude, shelf/coast geometry, and connected open basin.
6. Compute SST anomaly and heat transport:
   - conserve global mean heat anomaly.
   - smooth within basins, spread limited influence to adjacent coasts.
   - expose `ocean.sst_anomaly` for C4a/C4c/C4d.

Tests:

- Currents are zero outside solved ocean and tangent everywhere.
- Streamfunction contours do not cross land.
- Basins with no meaningful wind curl have weak gyres.
- Open subtropical basins produce coherent gyres rather than isolated coast
  patches.
- Straits limit exchange between basins and do not behave like open ocean.
- Upwelling occurs mainly on physically consistent eastern-boundary cold-current
  coasts.
- Heat transport has positive and negative regions but near-zero global mean.

Assets:

- `gyres.png`
- `current_streamfunction.png`
- `currents.png`
- `boundary_current_type.png`
- `ocean_heat_transport.png`
- `sst_anomaly.png`
- `upwelling.png`
- `strait_exchange.png`

Expected visual result:

- Ocean current maps read as basin-scale gyres and boundary currents rather than
  latitude bands or isolated coastline marks.

Implementation notes:

- `ClimateModule._ocean_currents` now solves a reduced graph-relaxation basin
  streamfunction from annual wind-stress curl plus planetary gyre structure,
  then derives tangent currents from streamfunction gradients.
- Existing compatibility fields remain: `ocean.currents`,
  `ocean.current_heat_transport`, `ocean.upwelling`, and `ocean.basin_id`.
- New C4b fields are emitted and archived:
  `ocean.current_streamfunction`, `ocean.gyre_id`,
  `ocean.boundary_current_type`, `ocean.strait_exchange`, and
  `ocean.sst_anomaly`.
- Validation now checks C4b field shape/finite status, streamfunction land
  leakage, gyre count, boundary-current area, strait exchange, and SST anomaly
  magnitude/mean.
- Current replay:
  `out_terminal_climate_replay_c4b1_20260705/`.
- Render probe:
  `out_terminal_climate_replay_c4b1_render_probe_20260705/earthlike_seed42/`
  writes `current_streamfunction.png`, `gyres.png`,
  `boundary_current_type.png`, `strait_exchange.png`, `sst_anomaly.png`, and
  `currents.png`.
- Six-world C4b metrics: ocean current speed p95 is about `0.17-0.21 m/s`;
  streamfunction ocean abs-p95 is about `0.99`; every world has nonzero gyre
  ids and SST anomaly structure.

### Phase C4c - SST, Wind, Evaporation, and Pressure Weak Coupling

Status: first pass complete; Earth comparison and all current Earth climate
gates pass on `out_terminal_climate_replay_c4c3_20260705/`

Purpose: make ocean currents and seasonal circulation exchange information
before precipitation is solved.  Warm/cold currents should affect SST,
evaporation, pressure, moisture source strength, and then weakly adjust winds.

Dependency note: C4c consumes C4a pressure/moisture diagnostics and C4b
streamfunction currents.  It should run as a bounded fixed-point iteration, not
as an unconstrained feedback loop.

Outputs:

- `climate.seasonal_sst`
- `ocean.sst_anomaly`
- `climate.ocean_heat_flux`
- Updated `climate.evaporation`
- Updated `atmosphere.seasonal_pressure_proxy`
- Updated `atmosphere.seasonal_wind`
- `climate.coupling_residual`

Algorithm:

1. Initialize seasonal SST from ocean cells in `climate.seasonal_temperature`.
2. Apply current heat transport and upwelling/cold-current cooling.
3. Convert SST to evaporation/source humidity with saturation-limited bulk
   formula and cold-current damping.
4. Feed warm/cold source-ocean pressure anomalies into seasonal pressure proxy.
5. Recompute bounded regional wind anomalies from the updated pressure field.
6. Recompute currents and SST anomaly for 2-4 iterations.
7. Stop when max bounded anomaly change is below threshold or iteration cap is
   reached.  Store `climate.coupling_residual`.

Tests:

- Coupling residual decreases or remains bounded across iterations.
- Global mean ocean heat anomaly remains close to zero.
- Warm-current regions have higher evaporation than comparable cold-current
  upwelling regions.
- Wind adjustment remains tangent and capped.
- Earthlike annual mean temperature remains within plausibility range.
- Frozen and waterworld presets remain numerically stable.

Assets:

- `seasonal_sst.png`
- `sst_anomaly.png`
- `ocean_heat_flux.png`
- `evaporation.png`
- `coupling_residual.png`

Expected visual result:

- Warm currents make local source oceans visibly more humid/active; cold eastern
  boundaries produce dry coasts before seasonal precipitation is added.

Implementation notes:

- `ClimateModule._weak_ocean_atmosphere_coupling` runs a bounded two-iteration
  loop over C4b SST anomaly/upwelling, ocean heat flux, pressure proxy,
  seasonal wind, and recomputed basin-streamfunction currents.
- `climate.seasonal_sst` and `climate.ocean_heat_flux` now feed seasonal
  evaporation, source-ocean warmth, seasonal pressure/moisture access, and
  final precipitation through the reduced hydroclimate path.
- Terminal replay arrays archive `climate.seasonal_sst`,
  `climate.ocean_heat_flux`, `climate.coupling_residual`,
  `climate.evaporation`, and `ocean.solved_mask`.
- Renderer assets now include `seasonal_sst.png`, `ocean_heat_flux.png`,
  `coupling_residual.png`, and `evaporation.png` alongside the C4b current/SST
  maps.
- Six-world C4c3 metrics: ocean heat flux remains mean-zero over the solved
  ocean mask, heat-flux abs-p95 is about `0.64-1.62 C`, coupling residual p95
  is `0.000-0.001 C`, and seasonal SST is clipped at seawater freezing
  (`-1.8 C`) over solved ocean cells.

### Phase C4d - Monsoon and Seasonal Hydroclimate

Status: first regionalization/object gate pass complete; C4d4 region-gate
replay preserves all current Earth climate gates

Purpose: replace static precipitation with four-season precipitation derived
from ITCZ/storm-track position, geography-derived pressure, terrain-steered
seasonal winds, SST-controlled source oceans, moisture access, and orographic
uplift/rain shadow.

Dependency note: C4d should not invent a separate monsoon rule.  Monsoon
rainfall should emerge where C4a/C4c diagnose summer heat lows, onshore flow,
warm source oceans, moisture corridors, and passable terrain.

Outputs:

- `climate.seasonal_precipitation`
- `climate.precipitation_seasonality`
- `climate.monsoon_index`
- `climate.dry_season_length`
- `climate.wet_season_peak`
- `climate.monsoon_rainfall_corridor`
- `climate.storm_track_rainfall_corridor`
- `climate.rain_shadow_index`
- `climate.regional_precipitation_response`
- `climate.hydroclimate_regions`
- Updated annual `climate.precipitation`
- Updated `climate.runoff`

Algorithm:

1. For each season, compute source humidity from seasonal SST, evaporation,
   sea-ice suppression, and upwelling/cold-current damping.
2. Advect/propagate moisture along seasonal wind and pressure-gradient flow,
   using moisture corridors from C4a.
3. Add convergence from:
   - ITCZ intensity and seasonal latitude.
   - storm-track intensity in winter hemisphere.
   - monsoon inflow where pressure gradient and moisture access align.
4. Add regional orographic precipitation:
   - use climate-scale barrier belts and windward normals.
   - distribute rain across windward zones, not one-cell ridge lines.
   - reduce leeward moisture through rain shadow corridors.
5. Aggregate four seasons:
   - annual `climate.precipitation` is the mean/equivalent total.
   - `climate.precipitation_seasonality` measures peak-season concentration.
   - `climate.monsoon_index` is signed summer-wet/winter-dry contrast on land.
   - `climate.runoff` uses seasonal or annual water balance depending on
     downstream readiness.

Tests:

- `climate.seasonal_precipitation` has shape `(4, n_cells)` and finite,
  nonnegative values.
- Annual precipitation equals the seasonal aggregate within tolerance.
- Large summer-heated continents with warm adjacent oceans and moisture access
  show positive monsoon index.
- Same continents have reduced or reversed monsoon signal in winter.
- Waterworld does not create fake land monsoons.
- Arid world keeps dry interiors but has plausible wet coasts or mountain
  windward zones where source water exists.
- Orographic precipitation is regional; high-relief concentration and
  precip-relief correlation remain below warning thresholds.
- Cold-current/upwelling coasts are drier than comparable warm-current coasts.
- Biomes and compiled hex maps still pass existing plausibility checks.

Assets:

- `precip_seasons.png`
- `precip_seasonality.png`
- `monsoon_index.png`
- `dry_season_length.png`
- `moisture_convergence.png`
- `orographic_precipitation.png`
- `monsoon_rainfall_corridor.png`
- `storm_track_rainfall_corridor.png`
- `rain_shadow_index.png`
- `regional_precipitation_response.png`
- `runoff.png`

Expected visual result:

- Tropical and subtropical rainfall follows moving ITCZ, monsoon corridors,
  source oceans, and terrain barriers rather than annual mountain stripes.
- Continents show wet coasts, dry interiors, and leeward deserts for causal
  geographic reasons.

Implementation notes:

- `ClimateModule._seasonal_hydroclimate` now derives four seasonal regional
  corridor fields from C4a/C4b/C4c diagnostics:
  `climate.monsoon_rainfall_corridor`,
  `climate.storm_track_rainfall_corridor`, `climate.rain_shadow_index`, and
  `climate.regional_precipitation_response`.
- Monsoon corridors consume monsoon potential, moisture access, thermal lows,
  onshore coast flow, warm source oceans, and terrain passability rather than
  using a separate rainfall shortcut.
- Storm-track corridors consume storm-track intensity, moisture access,
  onshore flow, and windward relief.  Rain-shadow corridors consume leeward
  exposure, barrier strength, wind gaps, and monsoon weakening.
- The regional precipitation response is applied season by season and then
  rescaled to preserve that season's land precipitation mean, so this pass
  changes spatial organization without acting as another broad rain boost.
- Terminal replay arrays archive the new fields plus `climate.runoff`, and the
  renderer writes the C4d corridor maps and `runoff.png`.

### Phase C5 - Ice, Snow, Cloud, and Vegetation Feedbacks

Status: pending

Purpose: soften hard climate boundaries and connect climate with biosphere.

Tasks:

- Make sea ice seasonal and persistent, then aggregate annual sea ice.
- Add simple snow persistence over land from seasonal temperature and precip.
- Add cloud/albedo proxy from humidity and storm-track intensity.
- Add optional vegetation albedo/evapotranspiration feedback from biome/NPP,
  with weak coupling to avoid runaway feedback loops.

Tests:

- Sea ice edge should not jump from 0 to 1 in one latitude band under Earthlike
  conditions.
- Frozen preset can still form broad ice cover.
- Vegetation feedback must not destabilize annual mean temperature.

Assets:

- Add `sea_ice_seasons.png`.
- Add `snow_ice.png`.
- Add `cloud_albedo.png` if useful.

Expected visual result:

- Polar and mountain climate boundaries become smoother and more physically
  interpretable.

### Phase C6 - Rendering, Compiler, and Archive Integration

Status: pending

Purpose: expose seasonal climate without breaking downstream systems.

Tasks:

- Keep annual fields as default for terrain erosion, biosphere, resources, and
  map compiler.
- Let biosphere optionally use seasonal minima/maxima:
  - frost limitation.
  - dry-season limitation.
  - growing season length.
- Let resources/weathering use annual runoff but keep seasonal runoff as a future
  diagnostic.
- Render climate sheets:
  - annual temperature/precip.
  - seasonal panels.
  - seasonality maps.
  - geography primitives and coupled climate diagnostics.
  - currents/winds.
- Extend `explain_cell` with climate seasonality story.

Tests:

- Existing compiled maps still run.
- Seasonal fields have correct shapes and finite values.
- `explain_cell` includes annual and seasonal climate fields when available.

Assets:

- `temperature.png`
- `temperature_seasons.png`
- `temperature_seasonality.png`
- `precip.png`
- `precip_seasons.png`
- `precip_seasonality.png`
- `geography_primitives.png`
- `ocean_basins.png`
- `terrain_barriers.png`
- `moisture_access.png`
- `monsoon_index.png`
- `sst_anomaly.png`
- `gyres.png`
- `currents.png`
- `wind_seasons.png`

### Phase C7 - Calibration Across Presets

Status: pending

Purpose: ensure climate improvements generalize beyond the Earthlike preset.

Presets to compare:

- `earthlike`: main qualitative comparison to Earth-style climate logic.
- `arid`: dry interiors, limited seas, strong continentality.
- `waterworld`: weak continentality, maritime climate, little monsoon.
- `frozen`: sea ice and snow-albedo feedback.
- `tidally_locked`: day-night circulation path should remain special.
- `stagnant_lid`: volcanic/high-CO2 world with different geography.

Acceptance checks:

- No preset crashes or produces NaNs.
- Earthlike remains within plausibility checks.
- Waterworld does not become fake monsoon land climate.
- Arid world keeps broad dry regions but not uniform desert everywhere.
- Frozen world keeps cold climate without single-band numerical walls.

## Reality Comparison Checklist

Use this checklist when visually comparing generated maps to Earth-like climate:

- Temperature is broadly latitude-controlled but not purely zonal.
- Ocean basins, shelves, straits, coasts, and terrain barriers are visible as
  reusable diagnostics rather than implied separately in each model.
- Same-latitude coastlines can differ because of currents and winds.
- Large land interiors show stronger seasonal amplitude than coasts.
- Monsoon-capable regions align with large heated continents, warm source
  oceans, onshore seasonal pressure gradients, and passable terrain corridors.
- Tropical rainfall follows ITCZ/monsoon patterns, not just mountains.
- Subtropical dry belts exist but are broken by coastlines, mountains, and monsoon.
- Mid-latitude storm tracks create wet belts on appropriate coasts.
- Orographic rain is regional and has leeward dry zones, not narrow cell stripes.
- Sea ice margins and snowlines are not hard one-cell walls.
- Biomes reflect seasonal stress, not only annual mean temperature and rain.

## Validation Plan

Add or update tests in `tests/test_engine.py` and diagnostics in
`aevum/validation.py`.

Candidate numeric gates:

- `max_adjacent_lat_band_delta_C < 18` for Earthlike annual mean.
- `land_seasonal_temp_amplitude_p50 > ocean_seasonal_temp_amplitude_p50`.
- `NH_land_JJA_mean_C > NH_land_DJF_mean_C` for Earthlike.
- `continent_components_cover_land == true`.
- `ocean_basin_components_cover_ocean == true`.
- `coast_orientation_max_normal_component < 1e-6`.
- `shelf_index_near_coast_p75 > shelf_index_deep_ocean_p75`.
- `strait_index_land_max == 0`.
- `barrier_index_high_relief_p75 > barrier_index_lowland_p75`.
- `moisture_access_leeward_mean < moisture_access_windward_mean` for strong
  barrier belts.
- `monsoon_potential_waterworld_p95 < small_threshold`.
- `monsoon_potential_summer_continent > monsoon_potential_winter_continent`
  where large warm continents have adjacent source oceans.
- `current_streamfunction_land_abs_max == 0`.
- `current_heat_transport_global_area_mean ~= 0`.
- `upwelling_eastern_boundary_fraction > upwelling_western_boundary_fraction`
  for open subtropical basins where geometry supports this comparison.
- `coupling_residual_p95` remains bounded after C4c iterations.
- `precip_seasonality_p75 > 1.2` on monsoon-capable Earthlike land.
- `currents_cross_land_cells == 0`.
- `coastal_temperature_asymmetry_index > small_threshold` where enough coasts
  exist.
- `finite_seasonal_fields == true`.
- `annual_precip_matches_seasonal_aggregate == true`.
- `cold_current_coast_precip_mean < warm_current_coast_precip_mean` after C4d
  where comparable coasts exist.
- `orographic_precip_concentration` and `precip_relief_correlation` stay below
  warning thresholds after seasonal precipitation is introduced.

Keep early thresholds permissive.  Tighten them only after assets improve.

## Progress Tracker

Legend:

- `[ ]` not started
- `[~]` in progress
- `[x]` complete
- `[!]` blocked or needs redesign

Milestones:

- [x] C0 diagnostics for current climate artifacts.
- [x] C1 seasonal insolation and land-ocean thermal inertia.
- [x] C2 seasonal winds, ITCZ migration, and storm tracks.
- [x] C2.5 geography-driven circulation anomalies.
- [x] C3 basin-aware ocean current heat transport.
- [x] C3.5 shared geography primitive layer.
- [x] C4a geography-derived seasonal pressure and moisture access.
- [x] C4b basin-streamfunction ocean currents.
- [x] C4c SST/wind/evaporation/pressure weak coupling.
- [~] C4d monsoon and seasonal hydroclimate.
- [!] C5 ice/snow/cloud/vegetation feedbacks blocked until Replay-R2-R6 pass.
- [~] C6 rendering/compiler/archive integration.
- [ ] C7 preset calibration and acceptance review.

Open issues:

- C4a consumes shared geography primitives for pressure, moisture access, and
  monsoon potential; C4b emits streamfunction currents plus SST anomaly; C4c
  now weakly couples SST, evaporation, pressure, and winds.  These pieces are
  still reduced proxies and need future calibration against monthly current/SST
  products, but the architectural coupling path is in place.
- Seasonal precipitation now has an initial C4d implementation, first-pass
  monsoon/storm/rain-shadow corridor fields, and multiple Earth gates.  The
  remaining C4d issue is upgrading these corridor proxies into more coherent
  regional object/flow structures with stronger visual organization.
- Current precipitation may still have orographic snake-band artifacts in some
  seeds; C4d diagnostics should quantify this against real Earth rain-shadow
  and coastal wetness envelopes.
- Current renderer supports seasonal precipitation, C4a/C4c/C4d diagnostic
  panels, runoff, and regional corridor maps, but the generated-vs-Earth
  contact sheet is still pending.
- C4b streamfunction currents are implemented, but the western/eastern boundary
  classification remains a reduced proxy; future OSCAR/monthly-current
  calibration can tighten it.
- Monsoon potential is now diagnosed from pressure, warm source oceans, terrain
  corridors, moisture access, and weak-coupled SST/current feedback; monsoon
  rainfall has first-pass C4d corridor regionalization but is not final until
  those corridors are organized as larger regional flow objects.
- Climate can be stale relative to final terrain in some extreme presets when
  terrain/ocean-mask drift does not cross the current re-solve trigger.

## Development Log

Use this section as an append-only log.  Include date, changed files, generated
assets, validation results, and remaining concerns.

Template:

```text
YYYY-MM-DD - Phase C?
Changed:
- ...
Validation:
- ...
Assets:
- ...
Concerns:
- ...
Next:
- ...
```

2026-06-21 - Plan archived
Changed:
- Added this climate system plan.
Validation:
- Not applicable; documentation only.
Assets:
- None.
Concerns:
- The current climate implementation remains annual-mean and lacks seasons,
  ocean heat transport, and true monsoon behavior.
Next:
- Start with C0 diagnostics, then C1 seasonal insolation and thermal inertia.

2026-06-21 - Phase C0 diagnostics implemented
Changed:
- Added `validation.climate_diagnostics()` and
  `validation.check_climate_diagnostics()`.
- Added diagnostics for adjacent latitude-band temperature jumps, neighbour-cell
  temperature jumps, land/ocean temperature contrast, orographic precipitation
  concentration, coastal temperature asymmetry, and seasonal-field availability.
- Added optional climate diagnostic rendering for scalar fields such as
  `climate.temperature_seasonality`, `climate.precipitation_seasonality`,
  `climate.monsoon_index`, and `ocean.current_heat_transport`.
- Added tests covering the current generated world and synthetic climate
  artifacts.
Validation:
- `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 24 tests in 140.36s.
Assets:
- No new assets generated by C0; rendering support is available for later fields.
Concerns:
- The warnings are expected to report missing seasonal climate fields and absent
  ocean heat transport until C1-C4 are implemented.
Next:
- Start C1 seasonal insolation and land-ocean thermal inertia.

2026-06-21 - Phase C1 seasonal temperature implemented
Changed:
- Added four-season `climate.seasonal_temperature` from seasonal insolation,
  eccentricity modulation, land-ocean thermal inertia, maritime lag, and
  continentality.
- Added `climate.temperature_seasonality` and `climate.continentality`.
- Kept `climate.surface_temperature` as the four-season mean for downstream
  compatibility.
- Added Feature Registry contracts and `temperature_seasons.png` rendering.
- Added tests for Earthlike seasonal asymmetry, land/ocean thermal inertia, and
  tidally locked non-seasonal behavior.
Validation:
- Targeted C1 and registry tests passed.
- `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 26 tests in 132.87s.
Assets:
- `temperature_seasons.png` and `temperature_seasonality.png` are now supported
  by rendering when a world is rendered.
Concerns:
- Seasonal precipitation, wind migration, monsoon behavior, and ocean heat
  transport remain future phases.
Next:
- Start C2 seasonal winds, ITCZ migration, and storm tracks.

2026-06-21 - Phase C1 annual EBM calibration
Changed:
- Rebased climate lapse cooling on elevation relative to sea level.
- Removed duplicate elevation cooling from the OLR solve; lapse cooling now
  enters once when converting the sea-level EBM state to surface temperature.
- Retuned the EBM OLR intercept from `210` to `195 W m^-2` to bring the
  Earthlike annual mean out of the severe cold bias.
- Reduced the seafloor-weathering scale from `0.65` to `0.45` so the carbon
  box does not overdraw Earthlike CO2 during warm calibration runs.
- Added Earthlike annual-mean temperature lower-bound checks.
Validation:
- Earthlike 2500-cell final sample: annual global mean `13.34 C`, land
  `14.78 C`, ocean `12.75 C`, seasonal min/max `-38.19/28.93 C`, CO2
  `156 ppm`.
- `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 26 tests in 127.37s.
Assets:
- No canonical assets regenerated in this step.
Concerns:
- CO2 remains lower than modern Earth, so the carbon box should be recalibrated
  later against degassing/weathering equilibrium rather than only temperature.
- Tropical and subtropical warm extremes remain muted until ocean heat transport,
  seasonal winds, and precipitation/cloud proxies are improved.
Next:
- Implement C2 seasonal winds, migrating ITCZ, and storm tracks.
- Then implement C3 ocean-current heat transport, because coastal temperatures
  and warm/cold-current asymmetry still cannot be fixed by EBM calibration alone.

2026-06-21 - Phase C2 seasonal circulation implemented
Changed:
- Added `atmosphere.seasonal_wind`, `atmosphere.itcz_latitude`,
  `atmosphere.itcz_intensity`, and `atmosphere.storm_track_intensity`.
- Replaced fixed annual wind bands with seasonally migrating Hadley/Ferrel/polar
  bands; annual `atmosphere.wind` remains the seasonal mean.
- Added winter-hemisphere storm-track strengthening and ITCZ migration tied to
  obliquity.
- Fixed the tidally locked wind projection so wind vectors remain tangent to the
  local sphere surface.
- Added Feature Registry contracts, climate diagnostics, tests, and rendering
  support for `wind_seasons.png` and `itcz_track.png`.
Validation:
- Targeted C2 and registry tests passed.
- Earthlike 1400-cell sample: ITCZ latitudes `[-15.939, 0.0, 15.939, 0.0]`,
  storm-track winter/summer ratios `1.84x` in both hemispheres, max wind normal
  component `2.7e-15`.
- `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 27 tests in 133.12s.
Assets:
- `wind_seasons.png` and `itcz_track.png` are now supported by rendering.
Concerns:
- Precipitation is still annual mean; C4 must consume the C2 fields to produce
  seasonal rainfall, monsoon, and dry-season behavior.
- Ocean currents still do not transport heat, so coastal asymmetry remains weak.
Next:
- Implement C2.5 geography-driven circulation anomalies.
- Then implement C3 ocean-current heat transport before C4, so seasonal
  hydroclimate can use both geography-aware winds and warmer/cooler source
  oceans.

2026-06-21 - Phase C2.5 planning inserted
Changed:
- Added Phase C2.5 between C2 and C3 for geography-driven circulation anomalies.
- Planned new fields for background wind, land-sea pressure proxy, thermal wind
  anomaly, orographic wind anomaly, and geographic circulation strength.
- Added algorithm steps for thermal pressure, coastal inflow/outflow, terrain
  steering, bounded wind composition, and downstream compatibility.
- Updated C3/C4 dependency notes so ocean currents and seasonal hydroclimate use
  geography-aware circulation rather than only zonal background winds.
Validation:
- Not applicable; documentation only.
Assets:
- Planned `thermal_wind_anomaly.png`, `land_sea_pressure.png`, and
  `geographic_circulation_index.png`.
Concerns:
- C2 remains a background circulation until C2.5 is implemented.
- C4 monsoon should not duplicate C2.5 logic; it should consume C2.5 fields.
Next:
- Implement C2.5 before C3/C4.

2026-06-21 - Phase C2.5 implemented and tested on new Earthlike layout
Changed:
- Implemented `atmosphere.background_seasonal_wind`,
  `atmosphere.land_sea_pressure_proxy`, `atmosphere.thermal_wind_anomaly`,
  `atmosphere.orographic_wind_anomaly`, and
  `atmosphere.geographic_circulation_index`.
- Composed final `atmosphere.seasonal_wind` from C2 background plus bounded
  thermal and orographic anomalies.
- Added Feature Registry contracts, climate diagnostics, rendering support, and
  tests for C2.5 fields.
- Tuned anomaly smoothing so thermal wind reads as regional geography response
  rather than a hard coastline outline.
Validation:
- Targeted C2.5 tests passed.
- New Earthlike layout `earthlike_c25_seed20260621`, seed `20260621`, 5000
  cells: land fraction `0.290`, annual mean `11.96 C`, mean precipitation
  `567 mm/yr`, CO2 `128 ppm`.
- Climate diagnostics passed; all validation checks passed for the generated
  world.
- C2.5 metrics on the generated world: background/final wind delta p95
  `2.14 m/s`, thermal anomaly p95 `1.95 m/s`, geographic circulation index p90
  `0.286`.
- `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 29 tests in 160.42s.
Assets:
- Generated at `/Users/rayw/Projects/aevum/out_c25_earthlike_seed20260621/`.
- Key new assets: `land_sea_pressure.png`, `thermal_wind_anomaly.png`,
  `orographic_wind_anomaly.png`, `geographic_circulation_index.png`,
  `wind_seasons.png`.
Assessment:
- C2.5 now makes wind fields respond to actual continents, coasts, and terrain.
- The new layout is useful for testing because it has large continents, open
  ocean basins, high-latitude land, and interior/coastal contrasts.
- The pressure proxy behaves physically: winter continents form cold highs,
  summer continents form heat lows.
- Remaining issue: precipitation remains annual and still shows line/terrain
  artifacts because C4 has not consumed the seasonal wind/pressure fields.
- Remaining issue: same-latitude coastal temperature asymmetry is still weak
  because C3 ocean heat transport is absent.
Next:
- Run final full suite after the asset iteration.
- Implement C3 basin-aware ocean current heat transport.

2026-06-21 - Phase C3 basin-aware ocean current heat transport implemented
Changed:
- Replaced wind-scaled diagnostic `ocean.currents` with basin-constrained
  spherical gyre proxies, western-boundary warm currents, eastern-boundary cold
  currents, and coastal upwelling.
- Added `ocean.current_heat_transport`, `ocean.upwelling`, `ocean.basin_id`,
  and `ocean.solved_mask`.
- Fed ocean-current heat transport into seasonal temperature and annual surface
  temperature, and added modest cold-current drying to evaporation and
  precipitation.
- Added Feature Registry contracts, climate diagnostics, rendering support, and
  regression tests for C3 fields.
Validation:
- Targeted C3 tests passed:
  `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation tests/test_engine.py::test_climate_diagnostics_cover_current_world_and_detect_artifacts -q`
  passed: 3 tests in 8.52s.
- Full suite passed:
  `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 31 tests in 166.23s.
- Earthlike 1400-cell sample: annual mean `12.43 C`, ocean heat transport
  absolute p95 `0.66 C`, land/coast heat influence absolute p95 `0.25 C`,
  coastal same-latitude asymmetry `0.70 C`, no solved-land currents, and climate
  diagnostics passed.
- Canonical C3 Earthlike seed `20260621`, 5000 cells: land fraction `0.290`,
  annual mean `10.96 C`, mean precipitation `561 mm/yr`, CO2 `162 ppm`, ocean
  heat transport absolute p95 `0.763 C`, coastal same-latitude asymmetry
  `0.660 C`, all validation checks passed.
Assets:
- `currents.png`, `ocean_heat_transport.png`, `upwelling.png`, and
  `ocean_basin_id.png` are now supported by rendering.
- Generated at `/Users/rayw/Projects/aevum/out_c3_earthlike_seed20260621/`.
Concerns:
- C3 is still a reduced proxy; it creates plausible first-order boundary-current
  contrasts but does not model thermohaline circulation, ENSO-like variability,
  or detailed basin streamfunctions.
- Waterworld/frozen short samples can report a stale climate ocean mask relative
  to final terrain; this is now a warning, not a C3 topology hard failure.
- Seasonal precipitation is still absent and remains the main climate warning.
Next:
- Follow the refined geography-coupled plan below: implement C3.5 shared
  geography primitives before C4a-C4d.

2026-06-21 - Geography-coupled climate plan refined
Changed:
- Updated Current State to reflect completed C1/C2/C2.5/C3 work instead of the
  original annual-only baseline.
- Added the Geography-Coupled Climate Architecture section.
- Added planned shared fields and objects for continents, ocean basins, shelves,
  straits, coast orientation, terrain barriers, wind gaps, moisture access,
  monsoon potential, basin streamfunctions, SST anomaly, and seasonal
  precipitation.
- Replaced the old single C4 plan with C3.5, C4a, C4b, C4c, and C4d:
  shared geography primitives; geography-derived seasonal pressure/moisture
  access; basin-streamfunction currents; SST/wind/evaporation/pressure weak
  coupling; and monsoon/seasonal hydroclimate.
- Expanded validation gates and progress milestones around geography topology,
  moisture corridors, streamfunction currents, coupled SST feedback, and
  seasonal precipitation aggregation.
Validation:
- Documentation-only update; tests not run.
Assets:
- No new assets generated.  Planned assets now include
  `geography_primitives.png`, `continent_id.png`, `ocean_basins.png`,
  `shelf_strait.png`, `coast_orientation.png`, `terrain_barriers.png`,
  `wind_gaps.png`, `seasonal_pressure.png`, `moisture_access.png`,
  `monsoon_potential.png`, `gyres.png`, `current_streamfunction.png`,
  `sst_anomaly.png`, `precip_seasons.png`, and `monsoon_index.png`.
Concerns:
- C3 remains a useful transition layer but should not become the final ocean
  model; C4b should replace its coast-rule currents with basin streamfunctions.
- C4d should not implement monsoon as a separate rainfall shortcut; it should
  consume C3.5/C4a/C4b/C4c diagnostics.
Next:
- Implement C3.5 shared geography primitive layer first.

2026-06-21 - Phase C3.5 shared geography primitives implemented
Changed:
- Added shared geography primitive fields:
  `climate.continent_id`, `climate.continent_interiority`,
  `climate.coast_orientation`, `climate.coast_distance`,
  `climate.coast_strength`, `climate.coast_facing_east`,
  `ocean.shelf_index`, `ocean.strait_index`, `terrain.barrier_index`, and
  `terrain.wind_gap_index`.
- Added object summaries:
  `climate.continents`, `ocean.basins`, `climate.coastline_segments`,
  `ocean.straits`, and `terrain.barrier_belts`.
- Added Feature Registry contracts, render outputs, climate diagnostics, and
  tests for C3.5.
- Iterated strait and wind-gap diagnostics after visual review: straits now
  filter isolated one-cell noise, and wind gaps are strong enough to inspect as
  low-pass/corridor candidates around terrain barriers.
Validation:
- Targeted C3.5 tests passed.
- Full suite passed:
  `/Users/rayw/Projects/aevum/.venv/bin/python -m pytest -q`
  passed: 33 tests in 208.36s.
- Canonical C3.5 Earthlike seed `20260621`, 5000 cells: land fraction `0.290`,
  annual mean `10.96 C`, mean precipitation `561 mm/yr`, CO2 `162 ppm`, all
  validation checks passed.
- C3.5 metrics on the generated world: `16` continent objects, `11` ocean-basin
  objects, `26` coastline segments, shelf/deep-ocean contrast `0.654`,
  strait ocean p95 `0.022`, `22` strait objects, terrain barrier high/low
  contrast `0.591`, wind-gap p95 `0.157`, and `17` barrier-belt objects.
Assets:
- Generated at `/Users/rayw/Projects/aevum/out_c35_earthlike_seed20260621/`.
- Key new assets: `geography_primitives.png`, `continent_id.png`,
  `ocean_basin_id.png`, `continent_interiority.png`, `coast_orientation.png`,
  `shelf_index.png`, `strait_index.png`, `terrain_barriers.png`,
  `wind_gaps.png`, and `c35_geography_contact_sheet.png`.
Assessment:
- The primitive layer now gives later phases shared, inspectable inputs for
  continents, basins, coasts, shelves, straits, barriers, and wind gaps.
- Shelf and terrain barrier diagnostics look coherent and align with the
  generated terrain.
- Strait detection is conservative after filtering; it is suitable as a C4b
  gateway candidate layer but not yet a true hydraulic exchange model.
- Wind gaps are interpretable as candidate low corridors but should be
  validated against actual seasonal flow in C4a/C4c.
Concerns:
- C3.5 does not yet change climate behavior; precipitation is still annual and
  still carries the expected C4 warning.
- Continent and basin fragmentation reflect the current generated geography.
  Future terrain/plate changes may need object-scale merging rules for tiny
  islands, marginal seas, and polar land bridges.
Next:
- Implement C4a: geography-derived seasonal pressure and moisture access using
  the C3.5 primitives.  This was completed later in the
  "C4a monsoon/moisture gate pass" entry below.

2026-07-05 - Initial terminal climate/biome and Earth calibration track
Changed:
- Added initial C4a/C4d climate fields:
  `climate.seasonal_precipitation`, `climate.precipitation_seasonality`,
  `climate.monsoon_index`, `climate.dry_season_length`,
  `climate.wet_season_peak`, `climate.moisture_convergence`,
  `climate.orographic_precipitation`, `atmosphere.seasonal_pressure_proxy`,
  `atmosphere.moisture_access`, `atmosphere.monsoon_potential`,
  `atmosphere.source_ocean_warmth`, and `atmosphere.terrain_blocking`.
- Updated the static biome classifier to consume seasonal temperature and
  seasonal precipitation stress instead of annual means only.
- Added terminal post-processing for accepted terrain worlds:
  `aevum.diagnostics.terminal_climate_biome` and CLI command
  `aevum terminal-climate-biome`.
- Added real-Earth calibration plan:
  `docs/EARTH_CLIMATE_REFERENCE_CALIBRATION_PLAN.md`.
- Added first executable Earth reference track:
  `aevum.diagnostics.earth_climate_reference` and CLI command
  `aevum earth-climate-reference`.
Validation:
- Targeted C4 climate tests passed earlier in this implementation pass; rerun
  after Earth-reference wiring is required before treating this as stable.
Assets:
- Terminal six-world climate/biome assets were generated at
  `out_terminal_climate_biomes_20260705/`.
- Earth reference CLI can now render same-grid ETOPO5, WorldClim land
  temperature/precipitation, NOAA PSL seasonal wind/pressure, and lightweight
  Koppen-Geiger major classes.
- Expanded Earth reference assets were generated at
  `out_earth_climate_reference_r2_20260705/`, including
  `earth_reference_8000cells_contact_sheet.png` and
  `earth_reference_24000cells_contact_sheet.png`.
- `earth_reference_*cells.npz` now includes:
  `earth.elevation_m`, `earth.land_mask`, `earth.monthly_temperature_C`,
  `earth.monthly_precip_mm`, `earth.seasonal_temperature_C`,
  `earth.seasonal_precip_mm_yr_equiv`, `earth.annual_temperature_C`,
  `earth.annual_precip_mm`, `earth.dry_month_count`,
  `earth.monthly_wind_u10_v10`, `earth.seasonal_wind_u10_v10`,
  `earth.monthly_slp_hPa`, `earth.seasonal_slp_hPa`,
  `earth.annual_slp_hPa`, `earth.seasonal_slp_anomaly_hPa`,
  `earth.koppen_class`,
  `earth.koppen_major_class`, and `earth.biome_class_proxy`.
Concerns:
- At this R4-reference point C4b and C4c were still pending; C4b later gained
  the first-pass streamfunction implementation recorded below, while C4c remains
  pending.
- OSCAR currents, GloH2O high-resolution Koppen, and true
  ESA/MODIS/RESOLVE biome/land-cover sources are registered in the manifest but
  not yet parsed into arrays.
- The current `earth.biome_class_proxy` is a Koppen-derived calibration proxy,
  not a final observed land-cover/biome product.
- The current terminal climate workflow rebuilds deterministic terrain from the
  accepted plate stack, then post-processes climate/biome without modifying the
  plate system.
Next:
- Implement an ocean-current reference path.
- Then compare generated temperature, precipitation, monsoon, and biome maps
  against Earth envelopes before tuning C4a/C4d parameters.

2026-07-05 - Earth reference R3 Koppen/biome upgrade
Changed:
- Upgraded Earth Koppen reference to GloH2O / Beck et al. V3, defaulting to the
  1991-2020 0.1 degree GeoTIFF source.
- Added RESOLVE Ecoregions 2017 polygon sampling to produce true terrestrial
  `earth.resolve_biome_class` and `earth.resolve_ecoregion_id` arrays.
- Kept the older Koppen-derived `earth.biome_class_proxy` as a coarse fallback
  and comparison layer, not as the primary observed biome reference.
Validation:
- `tests/test_earth_climate_reference.py` passes.
- R3 assets generated at `out_earth_climate_reference_r3_20260705/`.
- 24000-cell R3 output samples 27 Koppen fine classes, all 14 RESOLVE biomes,
  and 616 RESOLVE ecoregion ids; RESOLVE covers about 98.6% of sampled land
  area.
Current Earth-reference limits:
- Ocean currents and SST are still missing from real-Earth calibration.
- ETOPO5 should still be cross-checked against ETOPO 2022.
- Land-cover products such as MODIS/ESA/Copernicus are still useful, but they
  should be treated as land-cover cross-checks rather than biome replacements.
Next:
- Implement numeric surface-current reference, preferably OSCAR if accessible
  without auth friction, otherwise NOAA GODAS as a no-account fallback.
- Add SST reference before final climate/biome parameter tuning.

2026-07-05 - Earth reference R4 SST/current calibration anchors
Changed:
- Added NOAA OISST v2 1991-2020 monthly long-term mean SST and sea-ice
  references to the Earth calibration pipeline.
- Added NOAA/AOML annual near-surface drifter current climatology v3 as the
  first executable numerical ocean-current reference.
- Added `earth.monthly_sst_C`, `earth.seasonal_sst_C`, `earth.annual_sst_C`,
  `earth.monthly_sea_ice_concentration_pct`,
  `earth.seasonal_sea_ice_concentration_pct`,
  `earth.annual_sea_ice_concentration_pct`,
  `earth.surface_current_u_v`, `earth.annual_surface_current_speed_m_s`, and
  `earth.annual_surface_current_direction_deg` to the reference `.npz`.
Validation:
- `tests/test_earth_climate_reference.py` passes.
- R4 assets generated at `out_earth_climate_reference_r4_20260705/`.
- 24000-cell R4 output: ocean SST mean `18.3 C`, tropical SST mean `27.0 C`,
  SST seasonal-amplitude p90 `5.48 C`, annual sea-ice area above 15%
  concentration `8.9%` of valid ocean area, current-speed p90 `0.26 m/s`.
Current Earth-reference limits:
- Current reference is annual drifter climatology, not full seasonal OSCAR.
  It is enough to calibrate first-order observed current-speed envelopes and
  major current corridors, but not enough for seasonal current anomalies.
- The public ArcGIS OSCAR service exposes 2001-2020 monthly U/V LERC slices;
  using all 240 months remains a higher-cost upgrade.
Next:
- Compare the existing six terminal generated worlds against R4 Earth
  envelopes before tuning climate/biome parameters.
- Then begin climate-side tuning with temperature/SST influence and moisture
  access, keeping plate/terrain frozen.

2026-07-05 - Earth fitting gate hardening
Changed:
- Added `earth_climate_guardrails.csv` as the explicit per-check output from
  `aevum earth-climate-fit-report`.
- Added `--fail-on-guardrail` so CI-style sweeps can fail only on guardrail
  failures while allowing warning states to continue.
- Guardrail reports now count skipped checks separately from pass/fail checks.
Validation:
- Regenerated
  `out_earth_climate_fitting_f1_oceanfloor_gate_20260705/` with the hardened
  gate.  Current verdict remains `pass_with_warnings`: 0 failures, 1 warning,
  0 skipped checks.
Next:
- Continue Earth-first fitting with the high-rainfall tail as the main residual
  climate issue before deeper biome tuning.

2026-07-05 - F4 high-rainfall tail gate pass
Changed:
- Attributed the remaining Earth gate warning to a weak annual high-rainfall
  tail in warm, high-access, moderate-convergence earthlike regions, especially
  `earthlike_seed42`.
- Added a warm convective tail term in `ClimateModule._seasonal_hydroclimate`
  so p90/p95 rain increases without broadly lifting dry interiors.
Validation:
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_f4tail2_20260705/` and rendered PNG assets to
  `out_terminal_climate_replay_f4tail2_render_20260705/`.
- Re-ran Earth comparison at
  `out_earth_climate_comparison_f4tail2_20260705/`.
- Re-ran the hardened fitting gate at
  `out_earth_climate_fitting_f4tail2_gate_20260705/`.
- Current gate verdict: `pass`, with 0 failures, 0 warnings, and 0 skipped
  checks.
- Earthlike p90 precipitation ratios are now `0.517` and `0.638`; arid median
  precipitation and waterworld island guardrails remain in bounds.
Next:
- Keep plate/terrain frozen and move to a stricter Earth-pattern calibration
  pass before final biome threshold tuning: wet tropics, monsoon margins,
  windward mountains, dry subtropics, and cold high-latitude biome envelopes.

2026-07-05 - Earth pattern gate and first pattern pass
Changed:
- Added `aevum.diagnostics.earth_climate_pattern_gate` and CLI command
  `aevum earth-climate-pattern-gate`.
- The new gate compares broad Earth envelopes instead of exact map overlap:
  wet-tropics p90 and wet-area fraction, dry-subtropical fraction, monsoon wet
  area, mountain wet-tail fraction, high-latitude cold/ice-tundra fractions,
  and forest/tropical biome envelope.
- Strengthened geography-conditioned subtropical subsidence in F4 so earthlike
  dry belts appear without drying waterworld islands.
- Added a low-latitude warm convective tail term to support wet-tropics p90.
- Added smooth land-only polar cooling after the ocean SST floor so high-lat
  land can cool without reintroducing subfreezing ocean artifacts.
- Fixed biome precedence so cold-dry land is not overwritten as desert.
Validation:
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_pattern7_20260705/` and rendered PNG assets to
  `out_terminal_climate_replay_pattern7_render_20260705/`.
- Re-ran scalar comparison at
  `out_earth_climate_comparison_pattern7_20260705/`.
- Re-ran scalar fitting gate at
  `out_earth_climate_fitting_pattern7_gate_20260705/`; verdict `pass`, 0
  failures, 0 warnings.
- Ran Earth pattern gate at
  `out_earth_climate_pattern_gate_pattern7_20260705/`; verdict
  `pass_with_warnings`, 0 failures, 2 warnings.
- Pattern7 earthlike metrics: wet-tropics p90 Earth ratios `0.501` and
  `0.570`; dry-subtropical fractions `0.280` and `0.180`; high-latitude
  cold/ice-tundra fractions `0.696/0.458` and `0.915/0.819`.
Next:
- Treat climate-side Earth fitting as broadly ready for the next phase at this
  gate strictness.  The remaining work is F5-focused: compare generated biome
  semantics against Koppen/RESOLVE and tune forest/tropical classification
  without using more rainfall as a proxy.

2026-07-05 - F5 spatial biome organization gate pass
Changed:
- Added `aevum.diagnostics.earth_climate_spatial_biome_gate` and CLI command
  `aevum earth-climate-spatial-biome-gate`.
- The new gate checks latitude-level biome organization rather than total area
  alone: low-tropical forest/tropical cover, tropical low-latitude
  concentration, subtropical desert belts, cool-midlatitude forest/desert
  balance, and high-latitude tundra/ice versus desert.
- Initial attribution showed that coarse biome envelopes passed while
  cool-midlatitude forest belts remained too weak and cool/high-latitude land
  was overclassified as desert.
- Updated `BiosphereModule._biomes` so cool climates use lower forest and
  desert precipitation thresholds, and high-latitude cold-dry land is
  classified as tundra/ice before desert.
Validation:
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_f5spatial1_20260705/` and rendered PNG assets to
  `out_terminal_climate_replay_f5spatial1_render_20260705/`.
- Re-ran scalar comparison at
  `out_earth_climate_comparison_f5spatial1_20260705/`.
- Re-ran scalar fitting, Earth pattern, coarse biome, and spatial-biome gates:
  all verdicts are `pass`, with 0 failures and 0 warnings.
- Spatial-biome earthlike metrics: cool-midlat forest/desert fractions
  `0.328/0.119` and `0.395/0.198`; high-lat tundra+ice/desert fractions
  `0.877/0.003` and `0.980/0.000`.
Next:
- Continue Earth-first fitting at the next granularity: Koppen-like seasonal
  subtype organization, mountain ecological zonation, and wet/dry side
  asymmetry.

2026-07-05 - F5 seasonal subtype gate pass
Changed:
- Added `aevum.diagnostics.earth_climate_seasonal_subtype_gate` and CLI command
  `aevum earth-climate-seasonal-subtype-gate`.
- The new gate compares Earth and generated worlds using the same four-season
  dry-quarter definition.  It checks low-tropical seasonal-dry subtype area,
  low/mid-latitude dry-quarter coverage, broad precipitation seasonality
  amplitude, and high-latitude long-dry-season guardrails.
- Added a low-latitude precipitation-seasonality redistribution in
  `ClimateModule._seasonal_hydroclimate`.  It preserves annual precipitation
  per cell while sharpening existing wet/dry seasonal contrast.
- Retuned seasonal tropical biome semantics so wet Aw/Am-like tropical climates
  remain in the coarse tropical class; only drier, very seasonal tropical
  margins are converted to grassland.
Validation:
- Pre-fix seasonal subtype gate on `f5spatial1` failed only on
  `low_tropics_dry_quarter_ge2_fraction` for both earthlike seeds.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_f5seasonal2_20260705/` and rendered PNG assets
  to `out_terminal_climate_replay_f5seasonal2_render_20260705/`.
- Re-ran scalar fitting, Earth pattern, coarse biome, spatial-biome, and
  seasonal-subtype gates: all verdicts are `pass`, with 0 failures and 0
  warnings.
- Seasonal-subtype earthlike metrics: low-tropics dry-quarter >=2 fractions
  `0.098` and `0.096` versus Earth `0.138`; low/mid-latitude seasonality p75
  `3.063` and `3.146` versus Earth `2.270`.
Next:
- Continue Earth-first fitting with mountain ecological zonation and wet/dry
  side asymmetry before moving to broader climate-system architecture work.

2026-07-05 - F5 mountain zonation gate pass
Changed:
- Added `aevum.diagnostics.earth_climate_mountain_zonation_gate` and CLI
  command `aevum earth-climate-mountain-zonation-gate`.
- The new gate compares Earth and generated worlds using the same elevation,
  temperature, precipitation, and biome fields.  It checks high-mountain alpine
  ecology, high-mountain desert excess, midlatitude high-mountain desert
  excess, mountain cooling, and mountain wet-tail envelopes.
- Updated `BiosphereModule._biomes` so cool high-elevation land is treated as
  alpine/tundra before arid classification.  This is a biome semantic repair
  only; frozen terminal terrain and broad precipitation are unchanged.
Validation:
- Pre-fix mountain-zonation gate on `f5seasonal2` failed on high-mountain
  desert excess for both earthlike seeds and high-mountain alpine ecology for
  `earthlike_seed909`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_f5mountain1_20260705/` and rendered PNG assets
  to `out_terminal_climate_replay_f5mountain1_render_20260705/`.
- Re-ran scalar fitting, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, and mountain-zonation gates: all verdicts are `pass`, with
  0 failures and 0 warnings.
- Mountain-zonation earthlike metrics: high-mountain alpine ecology fractions
  `0.576` and `0.533`; high-mountain desert fractions `0.235` and `0.200`.
Next:
- Continue Earth-first fitting with wet/dry side asymmetry and true
  windward/leeward moisture organization before broader climate-system
  architecture work.

2026-07-05 - Windward/leeward precipitation gate pass
Changed:
- Added `aevum.diagnostics.earth_climate_windward_leeward_gate` and CLI command
  `aevum earth-climate-windward-leeward-gate`.
- The new gate compares Earth and generated worlds by deriving mountain-slope
  windward/leeward classes from seasonal winds and topographic gradients, then
  checking annual and seasonal precipitation contrast.
- Terminal climate array archives now save seasonal wind, background/thermal/
  orographic wind anomalies, terrain barrier and wind-gap fields, ITCZ/storm
  track diagnostics, and orographic precipitation so wind/terrain/precipitation
  coupling can be audited from replay artifacts.
- `ClimateModule` now computes slope-wind exposure from climate-scale
  topography and seasonal winds.  Orographic precipitation is redistributed
  from leeward to windward slopes with seasonal land-mean preservation, so this
  fix strengthens wet/dry side asymmetry without broad precipitation lifting.
Validation:
- Pre-fix windward/leeward gate on `f5windbase` failed for `earthlike_seed42`:
  annual windward/leeward precipitation ratio `1.068` versus Earth `2.081`,
  and seasonal median `1.166` versus Earth `2.063`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_f5wind3_20260705/` and rendered PNG assets to
  `out_terminal_climate_replay_f5wind3_render_20260705/`.
- Re-ran scalar fitting, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates: all
  verdicts are `pass`, with 0 failures and 0 warnings.
- Windward/leeward earthlike metrics: annual ratios `1.309` and `1.712`;
  seasonal median ratios `1.574` and `1.837`.
Next:
- Continue from Earth-fitting into the broader climate-system architecture:
  shared geography primitives, basin-aware currents, and monsoon/moisture
  corridors derived from terrain and ocean layout.

2026-07-05 - C4a monsoon/moisture gate pass
Changed:
- Added `aevum.diagnostics.earth_climate_monsoon_moisture_gate` and CLI command
  `aevum earth-climate-monsoon-moisture-gate`.
- Terminal climate arrays now also persist shared geography primitives useful
  for C4 work: continent id/interiority, coast distance/strength, basin id,
  shelf index, and strait index.
- `ClimateModule` now damps continent interiority and seasonal pressure by
  absolute landmass scale.  This keeps Earthlike/arid continents active while
  preventing waterworld islands from becoming artificial continental heat lows.
Validation:
- The new gate intentionally failed on the old `f5wind3` replay with `3`
  waterworld failures, all on `waterworld_seed707`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4a1_20260705/`.
- New C4a gate verdict:
  `out_earth_climate_monsoon_moisture_gate_c4a1_20260705/` -> `pass`, with
  `0` failures and `0` warnings.
- Existing Earth pattern, coarse biome, spatial-biome, seasonal-subtype,
  mountain-zonation, and windward/leeward gates all still pass on the C4a1
  replay.
- Tests:
  `41 passed in 5.20s` for the Earth climate gate subset, and
  `8 passed, 44 deselected in 167.63s` for the climate-related `test_engine`
  subset.
Next:
- Start C4b: replace the current basin-constrained current proxy with
  basin-streamfunction currents before doing C4c SST/wind/evaporation/pressure
  weak coupling.

2026-07-05 - C4b basin-streamfunction ocean currents first pass
Changed:
- Replaced the C3-style latitude/coast current core with a reduced
  basin-streamfunction solve inside `ClimateModule._ocean_currents`.
- Added C4b fields:
  `ocean.current_streamfunction`, `ocean.gyre_id`,
  `ocean.boundary_current_type`, `ocean.strait_exchange`, and
  `ocean.sst_anomaly`.
- Added feature-catalog entries, validation metrics, terminal array archival,
  and renderer output for the new C4b fields.
Validation:
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4b1_20260705/`.
- C4a monsoon/moisture, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates all pass on
  the C4b1 replay with `0` failures and `0` warnings.
- C4b six-world diagnostics: current speed p95 `0.17-0.21 m/s`,
  streamfunction abs-p95 about `0.99`, nonzero gyre ids in every world, and
  SST anomaly abs-p95 about `0.63-1.58 C`.
- Targeted tests:
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`
  and `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`
  passed; broader C4b/C4a subset passed with `6 passed in 78.14s`.
- Render probe for `earthlike_seed42` wrote the C4b assets under
  `out_terminal_climate_replay_c4b1_render_probe_20260705/earthlike_seed42/`.
Next:
- Implement C4c: weakly couple `ocean.sst_anomaly`, evaporation, pressure,
  moisture-source strength, and seasonal winds while preserving the C4a/C4b
  gate pass state.

2026-07-05 - C4c SST/wind/evaporation/pressure weak coupling first pass
Changed:
- Added bounded weak ocean-atmosphere coupling in `ClimateModule`: C4b SST
  anomaly and upwelling now produce mean-zero `climate.ocean_heat_flux`,
  seasonal SST, pressure adjustments, wind adjustments, recomputed currents,
  and `climate.coupling_residual`.
- `climate.seasonal_sst` and `climate.ocean_heat_flux` now feed seasonal
  evaporation, source-ocean warmth, seasonal pressure/moisture access, and
  precipitation through the reduced hydroclimate path.
- Added feature-catalog entries, validation metrics, terminal array archival,
  and renderer output for `climate.seasonal_sst`,
  `climate.ocean_heat_flux`, and `climate.coupling_residual`; terminal arrays
  also persist `climate.evaporation` and `ocean.solved_mask`.
- Renderer now writes `seasonal_sst.png`, `ocean_heat_flux.png`,
  `coupling_residual.png`, and `evaporation.png`.
Validation:
- Py-compiled the touched modules:
  `aevum/modules/climate.py`, `aevum/validation.py`, `aevum/render.py`,
  `aevum/features.py`, and
  `aevum/diagnostics/terminal_climate_biome.py`.
- Targeted tests:
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`,
  and
  `tests/test_engine.py::test_climate_diagnostics_cover_current_world_and_detect_artifacts`
  passed together with `3 passed in 70.70s`; the core C4c/C4b test also passed
  after the final archive/render edit.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4c3_20260705/`.
- C4a monsoon/moisture, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates all pass on
  the C4c3 replay with `0` failures and `0` warnings.
- Earth comparison on C4c3 reports `earthlike flagged: 0`; earthlike scores are
  `0.46` and `0.33`.
- Earth climate fit report on C4c3 reports guardrail verdict `pass` with `0`
  failures and `0` warnings.
- C4c3 six-world diagnostics: heat-flux ocean mean is `0.000000`, heat-flux
  abs-p95 is about `0.64-1.62 C`, coupling residual p95 is `0.000-0.001 C`,
  solved-ocean seasonal SST minimum is `-1.8 C`, and ocean evaporation means
  are about `28.7-39.6 mm/yr`.
- Render probe for `earthlike_seed42` wrote the C4c assets under
  `out_terminal_climate_replay_c4c3_render_probe_20260705/earthlike_seed42/`.
Next:
- Continue with C4d regional seasonal hydroclimate: make monsoon rain belts,
  storm-track precipitation, and rain shadows consume the C4a/C4b/C4c fields
  as coherent regional flow objects rather than mostly per-cell responses.

2026-07-05 - C4d regional seasonal hydroclimate first pass
Changed:
- Added first-pass C4d regional hydroclimate fields:
  `climate.monsoon_rainfall_corridor`,
  `climate.storm_track_rainfall_corridor`, `climate.rain_shadow_index`, and
  `climate.regional_precipitation_response`.
- The C4d response consumes existing C4a/C4b/C4c diagnostics: monsoon
  potential, moisture access, seasonal pressure, onshore coast flow, warm
  source oceans, storm-track intensity, windward/leeward exposure, terrain
  barriers, and wind gaps.
- Seasonal precipitation is shaped by these regional wet/dry corridors and then
  rescaled to preserve each season's land precipitation mean, so the pass
  improves spatial organization without a broad precipitation lift.
- Added feature-catalog entries, validation metrics, terminal array archival,
  renderer panels, and test coverage for the new C4d fields.  Terminal replay
  arrays also persist `climate.runoff`, and renderer now writes `runoff.png`.
Validation:
- Py-compiled the touched climate, validation, render, features, terminal replay
  archive, and test files.
- Targeted C4d test passed:
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  with `1 passed in 24.28s`.
- Broader climate subset passed:
  `test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `test_ocean_currents_are_basin_constrained_and_transport_heat`,
  `test_cold_boundary_currents_suppress_local_evaporation`, and
  `test_climate_diagnostics_cover_current_world_and_detect_artifacts` with
  `4 passed in 94.40s`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4d2_20260705/`.
- C4a monsoon/moisture, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates all pass on
  the C4d2 replay with `0` failures and `0` warnings.
- Earth comparison on C4d2 reports `earthlike flagged: 0`; earthlike scores are
  `0.46` and `0.33`.
- Earth climate fit report on C4d2 reports guardrail verdict `pass` with `0`
  failures and `0` warnings.
- C4d2 six-world diagnostics: seasonal precipitation still aggregates exactly
  to annual precipitation (`agg_delta=0.000000`), regional response p05/p95 is
  about `0.84-1.03`, earthlike monsoon-corridor land p90 is about
  `0.027-0.034`, storm-track corridor land p90 is about `0.141-0.165`, and
  rain-shadow land p90 is about `0.200-0.205`.
- Render probe for `earthlike_seed42` wrote the C4d assets under
  `out_terminal_climate_replay_c4d2_render_probe_20260705/earthlike_seed42/`,
  including `monsoon_rainfall_corridor.png`,
  `storm_track_rainfall_corridor.png`, `rain_shadow_index.png`,
  `regional_precipitation_response.png`, `precip_seasons.png`, and
  `runoff.png`.
Next:
- Continue C4d by adding visual/object-level gates for coherent monsoon rain
  belts, storm-track wet coasts, and leeward dry regions; then decide whether
  the current extracted regions need a true flow-network/seasonal-basin
  representation.

2026-07-05 - C4d seasonal hydroclimate object layer
Changed:
- Promoted the existing C4d corridor fields into
  `climate.hydroclimate_regions` diagnostic objects.  Each object records kind,
  season, threshold, area fraction, centroid, intensity, precipitation,
  moisture access, monsoon potential, regional response, and dominant continent
  id where available.
- Object kinds currently include `monsoon_rainfall_corridor`,
  `storm_track_rainfall_corridor`, `rain_shadow_region`,
  `wet_regional_precipitation_response`, and
  `dry_regional_precipitation_response`.
- The object layer is diagnostic-only: it does not alter seasonal precipitation,
  temperature, winds, currents, biomes, or runoff.
- Added the object feature contract, validation metrics, terminal replay summary
  counts, and test coverage.
Validation:
- Py-compiled the touched climate, feature, validation, terminal replay, and
  test files.
- Targeted C4d test passed:
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  with `1 passed in 25.16s`.
- Broader climate subset passed:
  `test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `test_ocean_currents_are_basin_constrained_and_transport_heat`,
  `test_cold_boundary_currents_suppress_local_evaporation`, and
  `test_climate_diagnostics_cover_current_world_and_detect_artifacts` with
  `4 passed in 95.53s`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4d3_objects_20260705/`; scalar replay metrics
  are unchanged from C4d2.
- Hydroclimate-region object counts in that replay:
  arid seeds `179` and `161`, earthlike seeds `163` and `190`, and waterworld
  seeds `49` and `52`.  Waterworld largest region fractions stay small
  (`0.0025` and `0.0084`), while earthlike/arid worlds expose larger regional
  structure.
- C4a monsoon/moisture, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates all pass on
  the C4d3 object replay with `0` failures and `0` warnings.
- Earth comparison on C4d3 objects reports `earthlike flagged: 0`; earthlike
  scores are `0.46` and `0.33`.
- Earth climate fit report on C4d3 objects reports guardrail verdict `pass`
  with `0` failures and `0` warnings.
Next:
- Add object-level/visual C4d gates that judge region coherence, seasonal
  migration, windward-leeward placement, and waterworld/arid false-positive
  limits instead of only checking field presence and broad Earth envelopes.

2026-07-05 - C4d hydroclimate-region object gate
Changed:
- Terminal climate replay now writes per-world `hydroclimate_regions.json`
  alongside `terminal_climate_arrays.npz`, preserving the full C4d region
  objects for diagnostics and manual audit.
- Added `aevum.diagnostics.earth_climate_hydro_region_gate` and CLI command
  `aevum earth-climate-hydro-region-gate`.
- The new gate checks object archive presence, Earthlike multi-kind and
  multi-season object coverage, coherent region-size proxies, and waterworld/
  arid false-positive limits for broad monsoon/wet-response objects.
Validation:
- Py-compiled the touched replay, CLI, diagnostic, and test files.
- New unit tests passed:
  `tests/test_earth_climate_hydro_region_gate.py` with `4 passed in 0.39s`.
- Replayed all six frozen terminal worlds to
  `out_terminal_climate_replay_c4d4_regiongate_20260705/`; scalar replay
  metrics are unchanged from C4d2/C4d3.
- Hydroclimate-region gate:
  `out_earth_climate_hydro_region_gate_c4d4_20260705/`, verdict `pass` with
  `0` failures and `0` warnings.
- Actual object audit on C4d4: arid seeds have `179` and `161` objects,
  earthlike seeds have `163` and `190`, and waterworld seeds have `49` and
  `52`; waterworld largest monsoon-region fractions are only `0.0011` and
  `0.0029`.
- C4a monsoon/moisture, Earth pattern, coarse biome, spatial-biome,
  seasonal-subtype, mountain-zonation, and windward/leeward gates all pass on
  the C4d4 replay with `0` failures and `0` warnings.
- Earth comparison on C4d4 reports `earthlike flagged: 0`; earthlike scores are
  `0.46` and `0.33`.
- Earth climate fit report on C4d4 reports guardrail verdict `pass` with `0`
  failures and `0` warnings.
Next:
- Continue C4d with visual/placement gates: region maps should be checked for
  seasonal migration, windward/leeward placement against mountain barriers,
  and storm-track wet-coast geometry.  After those gates exist, decide whether
  current region extraction is sufficient or needs explicit seasonal
  flow-network/river-basin style objects.

2026-07-05 - C4d hydroclimate-region placement proxies
Changed:
- Extended `earth-climate-hydro-region-gate` beyond object presence/coverage.
  It now reads both `hydroclimate_regions.json` and
  `terminal_climate_arrays.npz`.
- Added placement/seasonal metrics for:
  `monsoon_field_lat_shift_jja_minus_djf`,
  `monsoon_object_lat_shift_jja_minus_djf`,
  `storm_track_abs_lat_weighted_median`,
  `storm_track_coast_distance_weighted_median`, and
  `rain_shadow_dry_response_corr`.
- Earthlike checks now require summer-hemisphere monsoon migration, midlatitude
  storm-track placement, coastal/moisture-corridor storm-track attachment, and
  rain-shadow alignment with dry response.  Waterworld checks reject strong
  continent-like monsoon migration.
Validation:
- Py-compiled the hydro-region gate and tests.
- Updated hydro-region gate unit tests passed:
  `tests/test_earth_climate_hydro_region_gate.py` with `4 passed in 0.32s`.
- Placement-proxy gate on the C4d4 replay:
  `out_earth_climate_hydro_region_gate_c4d5_placement_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Key C4d5 placement metrics: earthlike monsoon field migration is `66.61` and
  `59.44` latitude degrees; earthlike monsoon object migration is `60.49` and
  `61.14`; earthlike storm-track weighted latitude is `39.52` and `38.31`;
  earthlike storm-track coast-distance medians are `0.058` and `0.078`; and
  rain-shadow/dry-response correlation is `0.978` and `0.983`.
- Waterworld false-positive guard remains clear: monsoon field migration is
  `17.86` and `23.91`, below the `35` threshold.
- Existing C4d4 Earth gates remain the scalar/pattern evidence for the same
  replay and still report `pass` with `0` failures and `0` warnings; no climate
  solver or biome logic changed in this C4d5 gate-only step.
Next:
- The remaining C4d validation gap is genuinely visual/map-based: render or
  sample object masks against seasonal precipitation panels to check that
  region shapes read as coherent belts/corridors rather than only passing
  scalar placement proxies.

2026-07-05 - C4d hydroclimate-region map-readability gate
Changed:
- Extended `earth-climate-hydro-region-gate` from placement proxies into a
  first map-readability gate.  The gate reconstructs graph adjacency from the
  archived `lat/lon` arrays and measures active-area fraction, largest connected
  component share, and boundary roughness for monsoon, storm-track, rain-shadow,
  wet-response, and dry-response maps.
- Earthlike checks now require monsoon and storm-track maps to be readable
  connected belts rather than isolated speckles.  Waterworld and arid checks
  reject broad false-positive monsoon or wet-response maps.
- Added a unit test that fails an earthlike run with a too-sparse monsoon map,
  so the new checks cover a real map-level failure mode instead of only object
  archive presence.
Validation:
- Py-compiled the hydro-region gate and tests.
- Updated hydro-region gate unit tests passed:
  `tests/test_earth_climate_hydro_region_gate.py` with `5 passed in 0.31s`.
- Map-readability gate on the C4d4 replay:
  `out_earth_climate_hydro_region_gate_c4d6_mapread_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Key C4d6 earthlike metrics: monsoon largest connected shares are `0.390` and
  `0.407`; monsoon active land fractions are `0.219` and `0.217`; storm-track
  largest connected shares are `0.344` and `0.190`; storm-track
  boundary-per-active-cell values are `1.891` and `1.487`; wet-response largest
  connected shares are `0.289` and `0.207`; dry-response largest connected
  shares are `0.181` and `0.314`.
Next:
- The remaining C4d validation gap is now rendered visual review/contact sheets
  and, if those reveal residual locality, explicit seasonal flow-network
  objects for moisture corridors and dry-season structure.

2026-07-05 - C4d hydroclimate-region contact sheets
Changed:
- `earth-climate-hydro-region-gate` now writes rendered contact sheets by
  default, with `--no-contact-sheet` available for fast metric-only runs.
- Each contact sheet is generated from the archived `terminal_climate_arrays.npz`
  used by the gate and shows four-season panels for seasonal precipitation,
  monsoon corridor mask, storm-track corridor mask, rain-shadow mask, and
  regional-response anomaly.
- The report JSON and Markdown now list emitted contact sheet paths.
Validation:
- Py-compiled the hydro-region gate, CLI, and tests.
- Updated hydro-region gate unit tests passed:
  `tests/test_earth_climate_hydro_region_gate.py` with `5 passed in 1.70s`.
- Contact-sheet gate on the C4d4 replay:
  `out_earth_climate_hydro_region_gate_c4d7_contacts_20260705/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Wrote `6` nonblank PNG contact sheets under
  `out_earth_climate_hydro_region_gate_c4d7_contacts_20260705/contact_sheets/`.
  A lightweight image audit confirmed all six PNGs are `1895 x 1380` and have
  varied pixels.
Next:
- C4d is considered sufficiently evidenced for the current Earth-fitting pass.
  The direct visual review found plausible seasonal migration and geographic
  coupling, with residual corridor-fragment behavior that should be handled in
  the next climate-system phase through explicit seasonal moisture-flow-network
  objects if we choose to deepen the hydroclimate architecture.

2026-07-05 - Earth-fitting pass accepted after C4d7 visual review
Changed:
- Reviewed the six C4d7 hydroclimate contact sheets directly after the metric
  gate.  Earthlike seeds show seasonally migrating monsoon wet cores,
  mid/high-latitude storm-track wet corridors, rain-shadow masks that align
  with dry regional response, and geography-coupled response anomalies.
- Arid and waterworld seeds do not show broad false-positive Earthlike monsoon
  systems; their wet corridors remain bounded or sparse, consistent with the
  C4d7 gate metrics.
- Promoted the residual issue from Earth-fitting blocker to climate-system
  design debt: monsoon and regional-response masks still read as corridor
  fragments rather than explicit routed moisture-flow-network objects.
Validation:
- Existing C4d7 hydroclimate-region gate evidence remains:
  `out_earth_climate_hydro_region_gate_c4d7_contacts_20260705/`, verdict
  `pass`, `0` failures, `0` warnings, `0` skipped checks, and `6` nonblank
  contact sheets.
- The Earth-fitting plan now records this acceptance so future work can use the
  gate suite as regression protection while moving back to climate-system
  development.
Next:
- Resume climate-system development after the accepted Earth-fitting pass.  The
  leading candidate is explicit seasonal moisture-flow-network objects; the
  alternative is to move to the next feedback layer only if flow-network
  structure is judged unnecessary for the next generated-world evaluation.

2026-07-05 - C4e seasonal moisture-flow-network first pass
Changed:
- Added a first explicit C4e moisture-flow-network layer in
  `aevum.modules.climate`.
- New seasonal fields:
  `atmosphere.moisture_flow_source`,
  `atmosphere.moisture_flow_pathway`, and
  `climate.moisture_flow_network_id`.
- New object set: `climate.moisture_flow_networks`, with
  `monsoon_moisture_flow_network`,
  `storm_track_moisture_flow_network`, and
  `mixed_moisture_flow_network` objects.
- The C4e extraction advects source-ocean moisture downwind through passable
  terrain, then intersects that pathway with C4d monsoon/storm/rain-shadow
  structure.  It is diagnostic/object-level only and preserves the accepted
  seasonal precipitation budget.
- Terminal climate archives now write the new C4e arrays plus
  `moisture_flow_networks.json`; the Feature Registry and climate validation
  diagnostics now expose the new fields and objects.
Validation:
- Py-compiled the touched climate, archive, feature, validation, and test files.
- Fast tests passed:
  `tests/test_climate_seasonal_redistribution.py` plus
  `tests/test_core.py::test_registry_resolves_dependencies` -> `3 passed`.
- Engine hydroclimate test passed:
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  -> `1 passed in 24.25s`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4e1_flow_20260705/`.
  All six worlds wrote `moisture_flow_networks.json` and C4e arrays.
- Existing C4d hydro-region regression gate on the C4e replay:
  `out_earth_climate_hydro_region_gate_c4e1_flow_20260705/`, verdict
  `pass`, with `0` failures, `0` warnings, and `0` skipped checks.
- Earth comparison and fitting guardrails on the C4e replay:
  `out_earth_climate_comparison_c4e1_flow_20260705/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4e1_flow_20260705/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C4e object sanity check: earthlike seeds produce four-season networks with
  maximum world-area fractions `0.026` and `0.0925`; waterworld maximum network
  areas remain small at `0.0064` and `0.0033`, so the new object layer does not
  create continent-scale false monsoon systems on waterworlds.
Next:
- C4e-specific moisture-flow-network gate/contact sheets are now implemented
  and archived below.  The next climate-system decision is whether C4e networks
  remain a diagnostic layer or become an active, conservative precipitation
  redistribution layer.

2026-07-05 - C4e moisture-flow-network gate and contact sheets
Changed:
- Added `aevum.diagnostics.earth_climate_moisture_flow_gate` and CLI command
  `aevum earth-climate-moisture-flow-gate`.
- The gate checks C4e object/archive presence, Earthlike multi-kind and
  four-season coverage, source-ocean moisture strength, routed land-pathway
  strength, pathway coupling to monsoon/storm support, pathway/network-id map
  readability, and waterworld/arid false-positive limits.
- The gate writes per-world C4e contact sheets with four-season panels for
  seasonal precipitation, source-ocean moisture, land moisture pathway,
  moisture-flow network id, and monsoon/storm support.
Validation:
- New unit tests passed:
  `tests/test_earth_climate_moisture_flow_gate.py` with `5 passed in 1.60s`.
- C4e gate on the six-world replay:
  `out_earth_climate_moisture_flow_gate_c4e2_gate_20260705/`, verdict `pass`
  with `0` failures, `0` warnings, `0` skipped checks, and `6` contact sheets.
- Image audit confirmed the six contact sheets are nonblank PNGs, all
  `1895 x 1380` with varied RGB channels.
- Direct representative visual review of `earthlike_seed42` and
  `waterworld_seed7` contact sheets found the intended structure: earthlike
  source moisture stays oceanic and enters land through seasonal pathways tied
  to monsoon/storm support, while waterworld pathways stay island-scale rather
  than expanding into broad continental monsoon networks.  Residual issue:
  flow-network id colors can still look segmented within broader corridors.
- Existing C4d hydro-region regression gate on the same C4e replay:
  `out_earth_climate_hydro_region_gate_c4e2_regression_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Targeted regression passed:
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  -> `14 passed in 25.02s`.
- Key C4e metrics on the two earthlike seeds: object counts `56` and `45`,
  four-season coverage, maximum network area fractions `0.026` and `0.0925`,
  land pathway p90 medians `0.949` and `0.953`, and pathway active-land
  fractions `0.335` and `0.359`.  Waterworld maximum network area fractions
  remain small at `0.0064` and `0.0033`, and pathway active-world fractions
  remain about `0.021` and `0.020`.
Next:
- Review the C4e contact sheets visually before any active precipitation
  coupling.  If accepted, prototype a conservative C4f precipitation
  redistribution pass where C4e moisture-flow networks modulate seasonal rain
  within each land/ocean budget instead of changing global totals.

2026-07-06 - C4f conservative moisture-flow precipitation response
Changed:
- Promoted C4e moisture-flow networks from diagnostic-only evidence into a
  conservative first-pass precipitation shaper in `aevum.modules.climate`.
- Added the seasonal field `climate.moisture_flow_precipitation_response`.
  Values below `1.0` mark land cells that donate seasonal precipitation, while
  values above `1.0` mark routed moisture-flow recipients.
- The C4f pass preserves each season's area-weighted land precipitation mean,
  leaves ocean precipitation unchanged, recomputes annual precipitation,
  evaporation, seasonality, monsoon index, dry-season length, and wet-season
  peak, then refreshes moisture-flow object mean precipitation from the shaped
  seasonal field.
- Terminal climate archives, feature metadata, validation diagnostics, and
  renderer output now include the C4f response.  Renderer response maps are
  centered on `1.0` so wetting and drying are visually separable.
Validation:
- Py-compiled the touched C4f files.
- Targeted tests passed:
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`.
- Render regression passed: `tests/test_render.py`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4f1_flowprecip_20260705/`.
- Earth comparison and fitting guardrails:
  `out_earth_climate_comparison_c4f1_flowprecip_20260705/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4f1_flowprecip_20260705/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing Earth climate gates all pass on the C4f replay with `0` failures,
  `0` warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a
  monsoon/moisture, C4d hydro-region, and C4e moisture-flow-network.
- Key C4f response metrics:
  earthlike seeds have land response p05/p95 around `0.807/1.068` and
  `0.794/1.069`; arid seeds around `0.805/1.072` and `0.751/1.074`; waterworld
  seeds are much weaker at about `0.880/1.040` and `0.895/1.033`, avoiding a
  broad false-positive continental monsoon response.  Maximum seasonal land
  mean deltas are numerical roundoff (`~1e-13 mm/yr`).
- Render probe:
  `out_terminal_climate_replay_c4f1_flowprecip_render_probe2_20260705/`.
  The `earthlike_seed42` response map is nonblank and centered on `1.0`.
Next:
- Keep plate and terrain frozen.  The next climate-system work should evaluate
  whether C4f needs basin/network-level moisture-budget constraints and better
  visual/object continuity before moving to stronger biome coupling.

2026-07-06 - C4f precipitation-response gate and contact sheets
Changed:
- Added `aevum.diagnostics.earth_climate_moisture_response_gate` and CLI command
  `aevum earth-climate-moisture-response-gate`.
- The gate checks C4f response archive presence, four-season shape, finite
  values, active diagnostic metadata, seasonal land-mean conservation,
  ocean-preserving response, bounded wet/dry strength, positive coupling to
  C4e moisture pathways and monsoon/storm support, rain-shadow avoidance,
  wet-response map readability, and waterworld false-positive limits.
- The gate writes metrics/check CSVs, a Markdown report, JSON summary, and
  contact sheets with seasonal precipitation, land moisture pathway, C4f
  response, flow-network id, monsoon/storm support, and rain-shadow panels.
Validation:
- New unit tests passed:
  `tests/test_earth_climate_moisture_response_gate.py` -> `5 passed`.
- C4f response gate without contact sheets:
  `out_earth_climate_moisture_response_gate_c4f2_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- C4f response gate with contact sheets:
  `out_earth_climate_moisture_response_gate_c4f2_contacts_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, `0` skipped checks, and `6`
  contact sheets.
- Image audit confirmed the six contact sheets are nonblank PNGs, all
  `1895 x 1601` with varied RGB channels.
- Direct visual review of representative `earthlike_seed42` and
  `waterworld_seed7` sheets found the intended pattern: earthlike response
  wetting/drying is tied to moisture pathways, monsoon/storm support, and
  rain-shadow structure, while waterworld response remains island-scale.
- Integrated C4d/C4e/C4f regression passed:
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  and `tests/test_render.py` -> `21 passed in 28.12s`.
Next:
- Keep using the accepted Earth gate suite as regression protection.  The next
  climate-system design step should decide whether to add basin/network-level
  conservative moisture budgets, because C4f currently preserves seasonal land
  means globally rather than by source basin or moisture-flow network.

2026-07-06 - C4g local moisture-budget precipitation conservation
Changed:
- Added the seasonal field `climate.moisture_budget_region_id`.
- C4f's active moisture-flow precipitation response now preserves seasonal
  precipitation means inside local land budget regions instead of preserving
  only the global land mean.  The first C4g authority is the climate continent
  component, so moisture-response wetting on one landmass cannot be paid for by
  drying an unrelated landmass.
- The implementation keeps plate and terrain frozen.  It changes only climate
  response shaping, field metadata, validation diagnostics, terminal replay
  archives, render output, and the C4f/C4g moisture-response gate.
- `earth-climate-moisture-response-gate` schema is now
  `aevum.earth_climate_moisture_response_gate.v2`.  The gate still checks the
  C4f response, and now also checks C4g budget-region archive presence, shape,
  finite values, region count, and per-budget mean conservation.
Validation:
- Targeted helper and gate tests passed:
  `tests/test_climate_seasonal_redistribution.py` and
  `tests/test_earth_climate_moisture_response_gate.py` -> `9 passed`.
- Integrated C4d/C4e/C4g regression passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 28.28s`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4g1_budget_20260706/`.
- C4g response gate:
  `out_earth_climate_moisture_response_gate_c4g1_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- C4g contact-sheet smoke:
  `out_earth_climate_moisture_response_gate_c4g1_smoke_contacts_20260706/`,
  verdict `pass` with `3` nonblank contact sheets.  Image audit sizes are
  `1895 x 1861` with varied RGB channels.
- Earth comparison and fitting guardrails on the C4g six-world replay:
  `out_earth_climate_comparison_c4g1_budget_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4g1_budget_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing Earth climate gates all pass on the C4g replay with `0` failures,
  `0` warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a
  monsoon/moisture, C4d hydro-region, C4e moisture-flow-network, and C4g
  moisture-response/local-budget.
- Key C4g smoke metrics:
  `earthlike_seed42` budget regions p50 `15`, max budget mean delta
  `4.55e-13 mm/yr`; `arid_seed101` p50 `1`, max delta `5.68e-14`; and
  `waterworld_seed7` p50 `22`, max delta `9.09e-13`.  Response p05/p95 remains
  bounded and waterworld response remains weak.
Next:
- Keep the accepted Earth fitting and C4g gates as regression protection.
- The next climate-system refinement can split the current continent-scale
  budget regions into explicit basin/network sectors where the source basin,
  pathway, and downwind receiver are stable enough to avoid over-damping the
  wet response.

2026-07-06 - C4h moisture-network sector budget refinement
Changed:
- `climate.moisture_budget_region_id` now starts from C4g continent-scale
  budget regions, then splits only large, coherent C4e moisture-flow networks
  into local C4h halo sectors.
- The split is deliberately conservative: small continents, small network cores,
  weak pathways, sectors that consume most of a continent, and tiny world-area
  fragments fall back to the parent continent budget.
- C4f/C4h response diagnostics now archive base budget-region count,
  post-split budget-region p50, sector split p50, and max per-budget seasonal
  precipitation mean delta.
- Terminal replay summaries now include C4h sector-split diagnostics, and the
  moisture-response gate schema is now
  `aevum.earth_climate_moisture_response_gate.v3`.
Validation:
- Targeted helper and gate tests passed:
  `tests/test_climate_seasonal_redistribution.py` and
  `tests/test_earth_climate_moisture_response_gate.py` -> `9 passed`.
- Integrated C4d/C4e/C4h regression passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 28.55s`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4h1_budget_sector_20260706/`.
- C4h response/local-budget gate:
  `out_earth_climate_moisture_response_gate_c4h1_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- C4h contact sheets:
  `out_earth_climate_moisture_response_gate_c4h1_contacts_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, `0` skipped checks, and `6`
  nonblank contact sheets.
- Earth comparison and fitting guardrails on the C4h six-world replay:
  `out_earth_climate_comparison_c4h1_budget_sector_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4h1_budget_sector_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing Earth climate gates all pass on the C4h replay with `0` failures,
  `0` warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a
  monsoon/moisture, C4d hydro-region, C4e moisture-flow-network, and C4h
  moisture-response/local-budget.
- Representative C4h split diagnostics from the replay: earthlike seeds have
  base budget regions `15` and sector-split p50 `7.5` / `3.5`; arid seeds have
  sector-split p50 `4.0` / `2.0`; waterworld seeds remain weak at `0.0` /
  `1.5`.  Max per-budget mean deltas remain roundoff-level
  (`~1e-13` to `9e-13 mm/yr`).
Next:
- Keep the accepted Earth fitting and C4h gate suite as regression protection.
- The next climate-system refinement should either tie budget sectors to
  explicit source-ocean basins/catchments, or move to precipitation object
  continuity if visual review shows the C4h sectors already read coherently.

2026-07-06 - C4i source-ocean basin attribution
Changed:
- Added `atmosphere.moisture_source_basin_id`, a four-season diagnostic field
  carrying the dominant ocean-basin id that supplies each moisture-flow pathway
  cell.  Ocean source cells inherit `ocean.basin_id`; land cells receive source
  labels through the same wind-aligned propagation used for C4e moisture
  pathways.
- Moisture-flow network objects now archive `source_basin_ids` and
  `dominant_source_basin_id`.
- C4h budget-sector splitting now respects source-basin attribution when it is
  available, so one local budget sector does not intentionally mix different
  diagnosed source oceans.
- Terminal replay archives, validation diagnostics, renderer output, and C4f
  moisture-response contact sheets now include the source-basin field.  The
  moisture-response gate schema is now
  `aevum.earth_climate_moisture_response_gate.v5`, and it explicitly checks
  source-basin archive presence, shape, finite values, attribution coverage for
  active pathways/wet response/networks, and source-basin purity inside budget
  sectors.
Validation:
- Targeted tests and registry/render checks passed:
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  and `tests/test_render.py` -> `12 passed in 27.13s`.
- Integrated C4d/C4e/C4i regression passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 29.04s`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4i1_source_basin_20260706/`.
- C4i moisture-response/local-budget/source-basin coherence gate:
  `out_earth_climate_moisture_response_gate_c4i2_source_basin_coherence_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- C4i contact sheets are nonblank; all six sheets are `1895 x 2095` and now
  include a source-ocean basin id row.
- Earth comparison and fitting guardrails on the C4i six-world replay:
  `out_earth_climate_comparison_c4i1_source_basin_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4i1_source_basin_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing Earth climate gates all pass on the C4i replay with `0` failures,
  `0` warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a
  monsoon/moisture, C4d hydro-region, C4e moisture-flow-network, and C4i
  moisture-response/local-budget.
- Representative source-attribution gate metrics from the replay: earthlike
  seeds attribute about `0.963` and `0.984` of seasonal land area to a source
  ocean basin; arid seeds attribute about `0.774` and `0.390`; active pathways,
  wet response corridors, and moisture-flow networks have source-basin
  attribution p50 `1.0`; budget source-purity p50 is `1.0` in all six worlds.
  Waterworld islands remain bounded by the existing waterworld false-positive
  response gate.
Next:
- Keep C4i as the current regression-protected state.
- The next climate-system step can move to precipitation object continuity.
  A stricter source-basin / receiver-catchment water-budget gate should wait
  until the precipitation objects expose stable receiver catchments; C4i2 now
  proves the current source-basin labels and budget sectors are coherent enough
  for that next object-continuity pass.

2026-07-06 - C4j precipitation-response object continuity
Changed:
- Added `climate.precipitation_response_region_id`, a four-season field that
  marks final wet/dry C4f precipitation-response objects after C4g/C4h/C4i
  local budget and source-basin constraints.
- Added `climate.precipitation_response_regions`, an object archive for final
  wet/dry response patches.  Each object records season, kind, area, centroid,
  response strength, precipitation, pathway, monsoon/storm support, rain
  shadow, dominant source basin, local budget region, and moisture-flow network
  attribution.
- C4j does not change the precipitation solve.  It is an object-continuity and
  auditability layer that prepares for later receiver-catchment budgets.
- Terminal replay archives, validation diagnostics, feature metadata, renderer
  output, and terminal summary JSON now include the C4j field/object layer.
- The moisture-response gate schema is now
  `aevum.earth_climate_moisture_response_gate.v6`; it checks C4j region-id
  archive presence/shape/finite values, object JSON presence, wet/dry kind
  coverage, source-basin attribution, budget-region attribution, and wet-object
  flow-network attribution.  Contact sheets now include a `C4j response region
  id` row.
Validation:
- Targeted object extraction tests passed:
  `tests/test_climate_seasonal_redistribution.py` -> `5 passed`.
- Gate tests passed:
  `tests/test_earth_climate_moisture_response_gate.py` -> `5 passed`.
- Integrated C4d/C4e/C4j regression passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `23 passed in 28.16s`.
- Six-world frozen terminal replay:
  `out_terminal_climate_replay_c4j1_precip_objects_20260706/`.
- C4j moisture-response/object-continuity gate with contact sheets:
  `out_earth_climate_moisture_response_gate_c4j2_contacts_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Earth comparison and fitting guardrails:
  `out_earth_climate_comparison_c4j1_precip_objects_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4j1_precip_objects_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing object regressions remain green on the same replay:
  `out_earth_climate_hydro_region_gate_c4j1_regression_20260706/` and
  `out_earth_climate_moisture_flow_gate_c4j1_regression_20260706/`, both
  verdict `pass`.
- Representative C4j gate metrics: earthlike seeds produce `71-83`
  precipitation-response objects across wet/dry kinds and four seasons, with
  p50 source-basin, budget-region, and wet-object flow-network attribution
  `1.0`; waterworld seeds produce `9-10` island-scale objects and remain within
  false-positive response bounds.
Next:
- Keep C4j as the current regression-protected state.
- The next climate-system step should reduce receiver-side fragmentation by
  merging C4j response objects into stable catchments/budget objects before
  enforcing a true source-basin -> receiver-catchment water-budget closure.

2026-07-06 - C4k receiver catchments and C5a2 circulation-layout reset
Changed:
- Added `climate.receiver_catchment_id` and `climate.receiver_catchments` as a
  diagnostic receiver-side catchment layer.  C4k merges C4j wet/dry response
  patches into seasonal land receiver catchments keyed by local budget region
  and source basin; it does not alter precipitation.
- Receiver catchments now include residual unbudgeted land/island cells as
  local residual accounting units, preventing small waterworld islands with
  `budget_id = -1` from disappearing from the receiver archive.
- Added `aevum.diagnostics.earth_climate_receiver_catchment_gate` and CLI
  command `earth-climate-receiver-catchment-gate`.  Waterworlds are no longer
  forced to contain broad C4j response-region bindings; their C4k checks focus
  on coverage, object archive, and island-scale bounds.
- After re-reading the original pipeline graph, the active work returned to the
  documented W/O order before further precipitation/biome tuning.  Added
  `aevum.diagnostics.earth_climate_circulation_layout_gate` and CLI command
  `earth-climate-circulation-layout-gate`.
- The new circulation-layout gate compares generated winds/currents against
  R6 Earth references and checks wind p90 magnitude, coastal onshore/offshore
  seasonal response, land-sea thermal wind anomaly, OSCAR-scale surface-current
  speed, and current land-mask leakage.
- `_seasonal_circulation` background winds were recalibrated to near-surface
  Earth wind magnitude.  The previous zonal template produced generated wind
  p90 about twice the R6/NCEP Earth seasonal 10 m wind p90, even though the
  downstream current-speed export was already OSCAR-scaled.
- Moisture-budget regions now have a source-basin fallback split for residual
  sectors, and arid gate semantics now distinguish active source-attributed
  pathways from huge weakly attributed dry interiors.
- Fixed a pre-existing one-line indentation blocker in `aevum/modules/terrain.py`
  that prevented later CLI imports; no terrain logic was changed.
Validation:
- Targeted tests passed:
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and targeted
  `tests/test_climate_seasonal_redistribution.py` cases -> `13 passed` plus
  source-budget tests -> `15 passed` in the focused runs.
- Seasonal hydroclimate/circulation regression passed:
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_earth_climate_monsoon_moisture_gate.py`,
  `tests/test_earth_climate_windward_leeward_gate.py`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  -> `9 passed in 38.86s`.
- C5a2 frozen six-world climate replay:
  `out_terminal_climate_replay_c5a2_circulation_sourcebudget_20260706/`.
- C5a2 circulation-layout gate:
  `out_earth_climate_circulation_layout_gate_c5a2_20260706/`, verdict `pass`
  with `0` failures, `0` warnings, and `0` skipped checks.
- C5a2 Earth comparison and fitting:
  `out_earth_climate_comparison_c5a2_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5a2_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C5a2 rendered replay/contact sheet:
  `out_terminal_climate_replay_c5a2_render_20260706/` and
  `out_earth_climate_comparison_c5a2_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- C5a2 downstream gates all pass with `0` failures and `0` warnings:
  `out_earth_climate_monsoon_moisture_gate_c5a2_20260706/`,
  `out_earth_climate_windward_leeward_gate_c5a2_20260706/`,
  `out_earth_climate_hydro_region_gate_c5a2_20260706/`,
  `out_earth_climate_moisture_flow_gate_c5a2b_20260706/`,
  `out_earth_climate_moisture_response_gate_c5a2b_20260706/`, and
  `out_earth_climate_receiver_catchment_gate_c5a2_20260706/`.
Next:
- Continue in the original dependency order.  The next step is not more
  precipitation or biome tuning; it is F2 spatial ocean-current/SST structure
  against OSCAR/OISST, followed by F3 pressure/wind/moisture placement.

2026-07-06 - C5b1 F2 ocean-current/SST spatial structure
Changed:
- Added `aevum.diagnostics.earth_climate_ocean_spatial_gate` and CLI command
  `earth-climate-ocean-spatial-gate`.
- The new gate uses R6 OSCAR/OISST fields to check current p90 magnitude,
  current land leakage, current-speed zonal dominance, SST zonal dominance,
  coastal swift-current SST anomaly spread, mean-zero heat transport, and
  earthlike near-coast/far-ocean strong-current placement.
- `_ocean_currents` now slightly strengthens diagnosed warm/cold boundary
  current vectors and mildly damps remote open-ocean current strength where
  there is little shelf, strait, or boundary-current influence.
- The change is deliberately limited to the ocean-current spatial structure
  layer.  Plate/terrain generation, temperature calibration, and precipitation
  coefficients are unchanged.
Validation:
- C5b1 six-world frozen replay:
  `out_terminal_climate_replay_c5b1_ocean_spatial_20260706/`.
- C5b1 rendered replay/contact sheet:
  `out_terminal_climate_replay_c5b1_render_20260706/` and
  `out_earth_climate_comparison_c5b1_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- F2 ocean-spatial gate:
  `out_earth_climate_ocean_spatial_gate_c5b1_20260706/`, verdict `pass` with
  `0` failures, `0` warnings, and `0` skipped checks.
- The original C5a2 issue was localized to `earthlike_seed42`: strongest
  current near-coast share `0.332` and far-ocean share `0.464`.  C5b1 moves
  those to `0.407` and `0.389`.  `earthlike_seed909` remains passing with
  near-coast share `0.520` and far-ocean share `0.266`.
- C5b1 circulation-layout gate:
  `out_earth_climate_circulation_layout_gate_c5b1_20260706/`, verdict `pass`.
- C5b1 Earth comparison/fitting:
  `out_earth_climate_comparison_c5b1_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5b1_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C5b1 downstream gates all pass with `0` failures and `0` warnings:
  monsoon/moisture, windward/leeward, hydro-region, moisture-flow,
  moisture-response, and receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  and `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`
  -> `8 passed in 53.01s`.
Next:
- Continue in the original dependency order with F3 pressure/wind/moisture
  placement.  Because temperature/SST, pressure/wind, currents, evaporation,
  moisture routing, and precipitation are mutually reinforcing, the next pass
  should add a coupled-consistency diagnostic before any isolated precipitation
  or biome tuning.  The current implementation is a bounded weak-coupling
  solver, not a full GCM, so the immediate target is internally consistent
  reduced physics rather than independent one-field matching.

2026-07-06 - C5c2 F3 coupled pressure/wind/moisture consistency
Changed:
- Added `aevum.diagnostics.earth_climate_coupled_consistency_gate` and CLI
  command `earth-climate-coupled-consistency-gate`.
- The gate checks the reduced model as a coupled system: seasonal warm land
  should correspond to low pressure proxy; winds should broadly follow
  high-to-low pressure gradients; SST should drive source-ocean warmth and
  evaporation; moisture/monsoon/storm/ITCZ support should explain seasonal
  precipitation; high monsoon-potential cells should also have adequate
  moisture and enhanced rainfall; cold-current/upwelling source regions should
  evaporate less than warm-current regions; precipitation and heat budgets
  should remain closed; and the weak ocean-atmosphere coupling residual should
  stay bounded.
- `_seasonal_pressure_moisture` now gates monsoon potential more strongly by
  moisture access.  This keeps dry thermal lows from becoming top monsoon
  potential unless they are also supplied by source-ocean moisture.
Validation:
- C5c2 six-world frozen replay:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`.
- C5c2 rendered replay/contact sheet:
  `out_terminal_climate_replay_c5c2_render_20260706/` and
  `out_earth_climate_comparison_c5c2_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- F3 coupled-consistency gate:
  `out_earth_climate_coupled_consistency_gate_c5c2b_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key metric movement: C5b1 `earthlike_seed42` top monsoon-potential moisture
  ratio was `0.504`; C5c2 improves it to `0.852`.  `earthlike_seed909` remains
  passing at `1.281`.
- F2/F3 regression gates pass:
  `out_earth_climate_ocean_spatial_gate_c5c2_20260706/`,
  `out_earth_climate_circulation_layout_gate_c5c2_20260706/`,
  `out_earth_climate_monsoon_moisture_gate_c5c2_20260706/`, and
  `out_earth_climate_windward_leeward_gate_c5c2_20260706/`.
- Earth comparison/fitting guardrails:
  `out_earth_climate_comparison_c5c2_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5c2_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Downstream object gates pass:
  `out_earth_climate_hydro_region_gate_c5c2_20260706/`,
  `out_earth_climate_moisture_flow_gate_c5c2_20260706/`,
  `out_earth_climate_moisture_response_gate_c5c2_20260706/`, and
  `out_earth_climate_receiver_catchment_gate_c5c2_20260706/`.
- Targeted tests passed:
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_earth_climate_monsoon_moisture_gate.py`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  and `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`
  -> `14 passed in 74.55s`.
Next:
- Continue in documented Earth-fitting order with F4 seasonal hydroclimate
  placement.  Use C5c2 coupled-consistency as a guardrail before any further
  precipitation redistribution or biome-threshold work.

2026-07-06 - C5d1 F4 seasonal hydroclimate placement gate
Changed:
- Added `aevum.diagnostics.earth_climate_seasonal_hydro_placement_gate` and
  CLI command `earth-climate-seasonal-hydro-placement-gate`.
- The gate evaluates F4 as an internal weak-coupling placement problem.  Wet
  seasonal land cells must be supported by moisture access, monsoon rainfall
  corridors, storm-track rainfall corridors, ITCZ intensity, regional
  precipitation response, and moisture-flow response.  Dry seasonal land cells
  must be explainable by low moisture, rain shadow, or regional dry response.
  Wet-season phase must usually align with the season of maximum support, and
  seasonal precipitation must aggregate exactly to annual precipitation.
- No climate formula changed in this pass.  The accepted C5c2 six-world replay
  already satisfies the new placement gate, so C5d1 is a diagnostic/guardrail
  hardening step rather than a retuning step.
Validation:
- Input replay:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`.
- F4 seasonal-hydro placement gate:
  `out_earth_climate_seasonal_hydro_placement_gate_c5d1_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key earthlike metrics: `earthlike_seed42` wet-support p25 ratio `2.301`,
  support/precipitation correlation `0.747`, wet-season peak-support match
  `0.922`, and rain-shadow precipitation ratio `0.556`; `earthlike_seed909`
  has `1.877`, `0.701`, `0.882`, and `0.358`, respectively.
- Targeted tests passed:
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_earth_climate_moisture_response_gate.py`
  -> `19 passed in 2.52s`.
Next:
- Keep C5d1 as the F4 guardrail and rerun it together with F2 ocean-spatial
  and F3 coupled-consistency gates after any change touching SST/current
  feedback, seasonal winds, moisture access, precipitation placement, or biome
  thresholds.  A future full coupling pass should add bounded iteration among
  temperature/SST, wind/pressure, ocean currents, evaporation, and
  precipitation; C5d1 only verifies that the current reduced fields agree.

2026-07-06 - C5d1 R6 full Earth-fitting acceptance bundle
Changed:
- Consolidated the accepted Earth-fitting state into one R6/C5d1 acceptance
  bundle before returning to broader climate-system implementation.  No climate
  formula changed in this step.
Validation:
- R6 comparison:
  `out_earth_climate_comparison_c5d1_acceptance_r6_20260706/`, `earthlike
  flagged: 0`.
- R6 fitting:
  `out_earth_climate_fitting_c5d1_acceptance_r6_20260706/`, guardrail verdict
  `pass` with `0` failures and `0` warnings.  Low-priority F3/F5
  `needs_tuning` labels remain future improvement handles, not acceptance
  blockers.
- Complete acceptance suite on
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/` passes with
  `0` failures, `0` warnings, and `0` skipped checks: pattern, biome,
  spatial-biome, seasonal-subtype, mountain-zonation, windward/leeward,
  monsoon/moisture, circulation-layout, ocean-spatial, coupled-consistency,
  seasonal-hydro placement, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Targeted acceptance tests passed:
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `27 passed in 3.23s`.
Next:
- Earth fitting is accepted at the current gate strictness.  Resume the climate
  system plan with bounded coupling iteration among temperature/SST,
  pressure/wind, currents, evaporation/moisture, and precipitation.  The C5e6
  bundle recorded below is now the active regression baseline.

2026-07-06 - C5e1 bounded precipitation-pressure/wind feedback
Changed:
- Added a conservative hydroclimate feedback pass inside `ClimateModule.step`.
  The model first computes preliminary seasonal hydroclimate, derives a small
  precipitation-pressure feedback from warm supported wet cores and dry
  low-moisture subsidence cells, caps/smooths the resulting pressure anomaly,
  derives a small tangent wind anomaly, and then recomputes final
  hydroclimate.  This is the first step toward the planned bounded coupling
  loop in which precipitation also feeds pressure/wind instead of only
  consuming them.
- Added archived fields:
  `atmosphere.precipitation_pressure_feedback`,
  `atmosphere.hydro_coupled_wind_anomaly`, and
  `climate.hydro_coupling_residual`.
- `terminal_climate_arrays.npz` now writes those fields for replay/gate
  diagnostics.
Validation:
- C5e1 six-world replay:
  `out_terminal_climate_replay_c5e1_hydro_feedback_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e1_hydro_feedback_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5e1_hydro_feedback_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Full C5e1 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement, hydro-region,
  moisture-flow, moisture-response, and receiver-catchment.
- Archived feedback bounds remain conservative: earthlike pressure-feedback
  abs-p95 is about `0.027-0.029`, wind-anomaly p95 is about `0.061-0.070 m/s`,
  and hydro-coupling residual p95 is about `0.022`; waterworld feedback is
  effectively zero.
- Targeted regression tests passed:
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `31 passed in 119.54s`.
Next:
- Expand the single hydro feedback pass into a small 2-4 iteration bounded
  coupling loop only if the current acceptance bundle remains stable.  C5e3
  below adds the formal coupling-convergence gate and becomes the active
  regression baseline.

2026-07-06 - C5e3 coupling-convergence gate and waterworld feedback fix
Changed:
- Added `aevum.diagnostics.earth_climate_coupling_convergence_gate` and CLI
  command `earth-climate-coupling-convergence-gate`.
- Added `atmosphere.land_sea_pressure_proxy` to `terminal_climate_arrays.npz`
  so the gate can compare precipitation-pressure feedback with the pre-feedback
  seasonal pressure field.
- The first C5e2 gate run correctly failed waterworlds: tiny island rainfall
  produced strong normalized hydro wind anomalies and pressure feedback.  C5e3
  fixes this in `ClimateModule._hydroclimate_pressure_wind_feedback` by scaling
  the hydro pressure/wind feedback by exposed land fraction.  Earthlike and
  arid continental cases remain effectively unchanged, while tiny-island worlds
  cannot create planet-scale hydro pressure systems.
Validation:
- C5e3 six-world replay:
  `out_terminal_climate_replay_c5e3_coupling_convergence_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e3_coupling_convergence_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e3_coupling_convergence_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e3_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key feedback metrics: earthlike pressure-feedback abs-p95 is about `0.035`,
  wind-anomaly p95 is `0.134-0.149 m/s`, and hydro residual p95 is about
  `0.033`; waterworld pressure-feedback abs-p95 is now only
  `0.00010-0.00026` and waterworld wind-anomaly p95 is about `0.003 m/s`.
- Full C5e3 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_coupling_convergence_gate.py`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `31 passed in 31.94s`.
Next:
- C5e3 is now the regression baseline for climate-system coupling work.  The
  next implementation step can expand the current single hydro feedback pass
  into a small bounded iteration loop, but it must keep the C5e3
  coupling-convergence and Earth acceptance gates green.

2026-07-06 - C5e4 bounded hydro-feedback iteration loop
Changed:
- Replaced the one-shot hydro feedback solve with a 3-pass bounded
  `_seasonal_hydroclimate_feedback_loop` inside `ClimateModule.step`.
- Each pass recomputes seasonal hydroclimate, derives conservative
  precipitation-pressure and tangent wind feedback, damps the update, and caps
  wind magnitudes before the final hydroclimate solve.
- Added archived field `climate.hydro_feedback_iteration_delta` plus
  `hydro_feedback_iteration_delta_p95` and `hydro_feedback_iteration_count`
  diagnostics.  The coupling-convergence gate now checks iteration-delta p95
  and max bounds, including stricter waterworld limits.
Validation:
- C5e4 six-world replay:
  `out_terminal_climate_replay_c5e4_hydro_iteration_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e4_hydro_iteration_r6_20260706/` reports
  `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e4_hydro_iteration_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- Coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e4_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key feedback metrics: earthlike wind-anomaly p95 is `0.136-0.150 m/s`,
  hydro residual p95 is about `0.033`, and iteration-delta p95 is
  `0.0000033-0.0000050`; waterworld wind-anomaly p95 remains about
  `0.003 m/s`, with effectively zero iteration delta.
- Full C5e4 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_coupling_convergence_gate.py`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `31 passed in 36.27s`.
Next:
- Treat C5e4 as the current climate-system regression baseline.  The next
  coupling work should add bounded response from winds/SST/evaporation back
  into ocean-current and moisture-source structure, while retaining the same
  convergence diagnostics and acceptance suite.

2026-07-06 - C5e5 evaporation-SST feedback in the ocean-atmosphere loop
Changed:
- Expanded `_weak_ocean_atmosphere_coupling` to 3 bounded iterations.
- Added `_ocean_evaporation_heat_feedback`, which derives a small SST heat-flux
  correction from seasonal ocean evaporation, cold-current/upwelling
  suppression, and ocean heat state.  The correction is smoothed over the
  solved ocean mask and forced to zero area-weighted ocean mean before it
  modifies `climate.ocean_heat_flux`.
- Added archived field `climate.ocean_evaporation_feedback`.
- Extended `earth-climate-coupling-convergence-gate` to require the new field
  and check feedback p95/max amplitude plus zero area-weighted ocean mean.
Validation:
- C5e5 six-world replay:
  `out_terminal_climate_replay_c5e5_evap_sst_feedback_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e5_evap_sst_feedback_r6_20260706/` reports
  `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e5_evap_sst_feedback_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- Coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e5_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key feedback metrics: ocean evaporation-feedback abs-p95 is about
  `0.17-0.21 C`, abs-max is `0.22-0.29 C`, area-weighted ocean mean is
  effectively `0`, ocean coupling residual p95 stays below `0.001 C`, and
  hydro iteration-delta p95 remains below `0.000007`.
- Full C5e5 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_coupling_convergence_gate.py`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`,
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `33 passed in 106.54s`.
Next:
- Treat C5e5 as the current climate-system regression baseline.  The next
  implementation step should add a bounded wind-stress/current response
  diagnostic and gate, so changes in seasonal wind do more than indirectly
  perturb currents through the existing streamfunction proxy.

2026-07-06 - C5e6 explicit wind-stress/current response
Changed:
- Added `ocean.wind_stress_current_response` as an explicit field emitted by
  `_ocean_currents`.
- The response is tangent, ocean-only, aligned with annual wind stress, capped
  as a small fraction of wind speed, and included in the current vector as a
  bounded surface response.  An initial open-ocean coefficient made
  `earthlike_seed42` fail the ocean-spatial gate because too many swift
  currents became far-ocean bands; the final coefficient keeps the response
  visible without breaking boundary-current placement.
- Added replay archive support for `ocean.wind_stress_current_response`.
- Extended `earth-climate-coupling-convergence-gate` to require and check the
  new field: ocean p50 presence, ocean p95 cap, land max zero, wind-stress
  alignment, and response-to-wind p95 ratio.
- Extended climate diagnostics/validation to report shape, non-finites,
  land leakage, normal component, and p50/p95 response speeds.
Validation:
- C5e6 six-world replay:
  `out_terminal_climate_replay_c5e6_wind_stress_current_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e6_wind_stress_current_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e6_wind_stress_current_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- Coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e6_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key feedback metrics: earthlike wind-stress response p50 is
  `0.042-0.047 m/s`, p95 is `0.067-0.070 m/s`, land max is `0`, alignment p50
  is `1.0`, and response-to-wind p95 ratio is about `0.012`.
- Ocean-spatial gate passes after the open-ocean response coefficient was
  reduced.  Earthlike current-speed p90 remains about `0.90` of Earth
  reference; swift-current far-ocean share is `0.405` for seed42 and `0.269`
  for seed909.
- Full C5e6 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile:
  `py_compile` passes for edited modules and tests.  Non-engine climate gate
  tests pass with `27 passed in 41.26s`.  The combined targeted pytest reached
  `4 passed` before being interrupted after `544.26s`, with the interruption
  inside an engine terrain/climate path; the standalone cold-boundary engine
  test also exceeded `5 min` in the current loaded environment and was
  interrupted.  This remains an incomplete engine-test run rather than an
  acceptance-gate failure.
Next:
- Treat C5e6 as the previous climate-system regression baseline for the
  following C5e7 work.  The next implementation step should improve
  source-ocean-basin and receiver-catchment moisture accounting so
  precipitation response can be checked against diagnosed seasonal ocean-basin
  supply, not only local pathway support.

2026-07-06 - C5e7 source-basin receiver-catchment accounting
Changed:
- Added `climate.source_basin_supply_index`, a bounded four-season diagnostic
  that combines source-ocean moisture strength, dominant source-basin labels,
  and routed landward pathway support.
- Added `climate.receiver_catchment_supply_balance`, a bounded four-season
  diagnostic that compares receiver-catchment precipitation magnitude with
  diagnosed source-basin supply support.  This is a consistency ledger, not a
  strict water-conservation solve.
- Receiver-catchment objects now archive
  `mean_source_basin_supply_index`,
  `source_basin_supply_attributed_fraction`,
  `source_basin_supply_mass_fraction`,
  `supply_supported_precipitation_fraction`, and
  `precipitation_supply_balance`.
- `terminal_climate_biome`, validation diagnostics, feature catalog, and
  seasonal renders now include the new source/receiver accounting fields.
- `earth-climate-receiver-catchment-gate` is upgraded to schema v2 and now
  requires both accounting fields plus nonzero source-backed wet-response and
  object-level balance metrics.
Validation:
- C5e7 six-world replay:
  `out_terminal_climate_replay_c5e7_source_receiver_accounting_20260706/`.
- Receiver-catchment v2 gate:
  `out_earth_climate_receiver_catchment_gate_c5e7_source_receiver_accounting_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Existing acceptance gates also pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, and moisture-response.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e7_visual_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e7_visual_r6_20260706/` reports guardrail
  verdict `pass` with `0` failures and `0` warnings.
- C5e7 visual comparison:
  `out_earth_climate_comparison_c5e7_visual_r6_20260706/earth_vs_generated_climate_contact_sheet.png`
  now includes generated temperature, precipitation, biome, and current-speed
  panels for all six replay worlds.  The comparison diagnostic renders missing
  generated preview PNGs directly from archived replay arrays, so this did not
  change climate values.
- `earth-climate-fit-report` is upgraded to
  `aevum.earth_climate_fitting_report.v2`; F3 scoring now uses the same
  low/mid-latitude summer moisture/monsoon diagnostics as the dedicated
  monsoon/moisture gate.  C5e7 final fit-report phase statuses are all
  low-priority `watch`: F1 `0.35`, F2 `0.10`, F3 `0.00`, F4 `0.00`, F5 `0.00`.
- Key F3 final-report metrics for earthlike seeds `42/909`: low/mid-latitude
  summer moisture p75 is `0.837/0.843`, summer monsoon-potential p90 is
  `0.358/0.667`, and summer-minus-winter monsoon-potential p75 is
  `0.259/0.483`.
- Key C5e7 receiver metrics: earthlike source-supply-attributed land p50 is
  `0.56/0.67`, wet-response source-supply p50 is `0.87/0.91`, and receiver
  supply-balance land p50 is `0.69/0.64` for seeds `42/909`.
- Tests/compile:
  `py_compile` passes for edited climate, diagnostics, validation, rendering,
  feature catalog, and tests.  Target tests pass:
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`,
  `tests/test_earth_climate_receiver_catchment_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`.
Next:
- Treat C5e7 as the previous climate-system regression baseline.  C5e8 below
  supersedes it by feeding source-basin receiver accounting back into a bounded
  regional precipitation solve.  The visual contact sheet also keeps the
  reduced-model latitude-banding residual visible for later F2/F3 mechanics
  work.

2026-07-06 - C5e8 receiver-supply precipitation feedback
Changed:
- Added `climate.receiver_supply_precipitation_feedback`, a four-season
  response field that records the second conservative precipitation
  redistribution pass from source-basin supply and receiver-catchment balance.
- The new pass runs after C4f and after an initial C5e7 source/receiver
  accounting solve.  It then recomputes hydroclimate objects, moisture-flow
  object precipitation summaries, C4j precipitation-response objects, C4k
  receiver catchments, and final C5e7 accounting from the adjusted hydro state.
- The feedback is deliberately small and bounded: it is land-only, keeps ocean
  feedback exactly `1.0`, clips land response conservatively, and preserves
  every seasonal `climate.moisture_budget_region_id` mean.
- Terminal replay arrays, validation diagnostics, feature metadata, renderer
  output, and receiver-catchment gate metrics now include the new feedback
  field.
- `earth-climate-receiver-catchment-gate` is upgraded to schema v3 and checks
  C5e8 archive presence, shape/finite values, bounded land response, and ocean
  neutrality.
Validation:
- C5e8 six-world replay:
  `out_terminal_climate_replay_c5e8_receiver_supply_feedback_20260706/`.
- C5e8 feedback diagnostics: earthlike seeds `42/909` have land response
  p05/p95 about `0.938/1.018` and `0.937/1.017`; waterworld seeds remain
  smaller at about `0.970/1.011` and `0.978/1.007`; max land and local-budget
  mean deltas are numerical noise.
- Receiver-catchment v3 gate:
  `out_earth_climate_receiver_catchment_gate_c5e8_20260706/`, verdict `pass`
  with `0` failures, `0` warnings, and `0` skipped checks.
- Existing acceptance gates also pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, and moisture-response.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e8_receiver_supply_feedback_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e8_receiver_supply_feedback_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.  Phase
  statuses remain low-priority `watch`: F1 `0.35`, F2 `0.10`, F3 `0.00`,
  F4 `0.00`, F5 `0.00`.
- Tests/compile:
  `py_compile` passes for edited climate, validation, feature catalog,
  rendering, diagnostics, and tests.  Target tests pass:
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_earth_climate_receiver_catchment_gate.py`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`.
Next:
- Treat C5e8 as the previous climate-system regression baseline.  C5e9 below
  supersedes it by addressing the reduced-model ocean/SST latitude-banding
  residual without broad precipitation lifting.

2026-07-06 - C5e9 ocean/SST structure anti-zonal pass
Changed:
- Strengthened the reduced basin gyre, boundary-current, and upwelling heat
  anomaly terms in `ClimateModule._ocean_currents`, while keeping ocean heat
  transport ocean-confined and area-mean neutral.
- Reduced over-smoothing of the ocean heat anomaly so basin and boundary-current
  structure remains visible at map scale.
- Exported ocean surface temperature now uses the coupled
  `climate.seasonal_sst` for ocean cells instead of the weaker pre-coupling
  heat-transport projection.  Land still uses the regular seasonal temperature
  path plus bounded current/coastal influence.
- `earth-climate-ocean-spatial-gate` now records same-latitude current/SST
  residual amplitudes and adds earthlike checks for minimum SST residual
  structure relative to Earth and minimum heat-transport p95.  This directly
  protects against a visually plausible but mostly zonal SST map.
- `validation.climate_diagnostics` and terminal replay summaries now expose
  `seasonal_sst_zonal_residual_abs_p95_C`, ocean heat-flux p95, and
  current-heat-transport p95.
Validation:
- C5e9 six-world replay:
  `out_terminal_climate_replay_c5e9_ocean_structure_20260706/`.
- Key earthlike metrics: seeds `42/909` have SST same-latitude residual p95
  about `3.07/3.39 C`, ocean heat-flux p95 about `1.28/1.25 C`, and
  current-heat-transport p95 about `1.08/1.06 C`.  Mean temperature remains
  stable at about `15.33/13.86 C`, and land precipitation p50 remains inside
  the accepted C5e8 envelope at about `297/402 mm/yr`.
- C5e9 ocean-spatial gate:
  `out_earth_climate_ocean_spatial_gate_c5e9_20260706/`, verdict `pass` with
  `0` failures, `0` warnings, and `0` skipped checks.  Earthlike SST residual
  ratios to Earth are about `0.64/0.71`, clearing the new `0.58` floor.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e9_ocean_structure_r6_20260706/` reports
  `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e9_ocean_structure_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.  F1-F5 phase
  statuses remain low-priority `watch`: F1 `0.35`, F2 `0.10`, F3 `0.00`,
  F4 `0.00`, F5 `0.00`.
- Existing acceptance gates also pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile:
  `py_compile` passes for edited climate, validation, terminal replay summary,
  ocean-spatial diagnostics, and tests.  Target tests pass:
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`, and
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`
  (`12 passed` for the combined target run).
Next:
- Treat C5e9 as the current terminal-world regression baseline only.  Do not
  move to ice/snow/cloud/vegetation feedbacks yet.  The real-Earth replay track
  is reset below: first repair and revalidate boundary conditions, pressure/
  wind, ocean currents, and SST/energy closure, then moisture and precipitation,
  and only then downstream cryosphere and biome semantics.

2026-07-06 - Replay-R geoscience subgraph order reset
Changed:
- Replaced the old "next C5 feedback" direction with an explicit Replay-R0-R9
  causal order in `docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md`.
- The new order is: R0 replay harness/invariants, R1 boundary forcing, R2
  pressure/wind, R3 ocean dynamics, R4 SST/energy closure, R5 moisture source
  and transport, R6 precipitation/hydroclimate objects, R7 cryosphere/cloud/
  vegetation feedback, R8 climate classes/biomes, and R9 generated-world
  promotion.
- Added the ordering rule that downstream maps may be rendered and scored only
  as observers until their upstream support layers pass.  In particular, sea
  ice, snow, clouds, vegetation feedback, Koppen, and biome tuning are blocked
  until R2-R6 are accepted on real-Earth replay.
- Marked the exploratory real-Earth replay outputs
  `out_real_earth_climate_replay_f3a_20260706/`,
  `out_real_earth_climate_replay_f3b_20260706/`,
  `out_real_earth_climate_replay_f3c_20260706/`,
  `out_real_earth_climate_replay_f4_20260706/`,
  `out_real_earth_climate_replay_f4b_20260706/`, and
  `out_real_earth_climate_replay_f4b_render_20260706/` as diagnostic-only,
  invalid promotion baselines because they mixed downstream interpretation with
  still-unsettled upstream wind/ocean/SST mechanics.
- Kept the real-Earth replay harness and C5e9 terminal-world regression bundle
  as useful evidence, but separated "terminal-world regression baseline" from
  "real-Earth replay promotion baseline".
Validation:
- No plate or terrain code was changed.
- `py_compile` passes for `aevum/modules/climate.py`,
  `aevum/modules/biosphere.py`,
  `aevum/diagnostics/real_earth_climate_replay.py`, and `aevum/cli.py`.
- Target replay harness test passes:
  `tests/test_real_earth_climate_replay.py` -> `1 passed in 0.83s`.
- Reset no-render real-Earth replay:
  `out_real_earth_climate_replay_replay_r2_r4_reset_20260706/`.
  Summary: surface-temperature MAE `2.94 C`, land-precipitation MAE
  `459.4 mm/yr`, annual wind p90 replay/Earth `6.72/6.72 m/s`,
  ocean-current p90 ratio `0.981`, and validation `PASS`.
- The reset replay still warns that annual temperature has a steep adjacent
  latitude-band jump.  It also shows weak SON pressure-pattern correlation
  (`0.032`) and large surface-temperature zonal residual p95 (`12.05 C`), so
  the next repair should remain in R2/R4 foundation mechanics.
Next:
- Use the reset replay as the new foundation baseline and continue in order.
  The first implementation target is R2/R4: pressure/wind seasonal phasing plus
  temperature/SST energy-wall behavior.  R3 ocean-current metrics are close
  enough to remain guarded while R2/R4 are repaired.  Precipitation, sea ice,
  and biome maps should be inspected only as observer outputs until those
  foundation layers pass.

2026-07-06 - Replay-R2/R4 phase4 foundation checkpoint, not promoted
Changed:
- Kept the wind-driving base pressure conservative, but added a smaller
  forward-cooling seasonal phase correction inside the final pressure/moisture
  layer.  This improves the weak SON pressure replay without letting the
  pressure correction directly drive the ocean-current solver.
- Changed land ice/snow albedo in the EBM from a single permanent-ice value for
  all cold land to a split between ordinary seasonal snow-covered land and
  high-latitude high-elevation permanent-ice settings.  This addresses the
  North high-latitude lowland overcooling found in real-Earth replay while
  keeping permanent ice-cap albedo high.
- Added Earth-aware adjacent latitude-band temperature-jump metrics to
  `real-earth-climate-replay` and its test.  This prevents the Antarctic/
  Southern Ocean natural gradient from being treated as the same kind of
  generic heat-wall warning as an unrealistic generated-world wall.
Validation:
- `py_compile` passes for the edited climate and replay diagnostics.
- `tests/test_real_earth_climate_replay.py` passes.
- Real-Earth phase4 replay:
  `out_real_earth_climate_replay_replay_r2_r4_phase4_20260706/`.
  Surface-temperature MAE improves to `2.62 C`, land-temperature MAE to
  `3.83 C`, ocean-SST MAE is `2.13 C`, and validation passes.  Seasonal
  pressure correlations are `0.617/0.601/0.682/0.205` for DJF/MAM/JJA/SON.
  The replay max adjacent latitude-band jump is `18.58 C`, below the Earth
  reference value `20.46 C` at the same Antarctic/Southern Ocean transition.
- Six-world phase4 replay:
  `out_terminal_climate_replay_replay_r2_r4_phase4_20260706/`.
  Earth comparison reports `earthlike flagged: 0`; fit report has `0`
  guardrail failures and `3` warnings.
- Blocking guardrail:
  `out_earth_climate_ocean_spatial_gate_replay_r2_r4_phase4_20260706/`
  fails with `8` failures.  Current-speed maps are too latitude-banded for both
  earthlike and waterworld guardrails, and earthlike SST same-latitude residual
  ratios are below the C5e9 floor.  A phase6 attempt to reduce the planetary
  curl template was rejected because it worsened the gate and was reverted.
Next:
- Do not proceed to R5 moisture, R6 precipitation, R7 cryosphere/cloud/
  vegetation, or R8 biome tuning.  The next implementation target is R3:
  restore basin/coast/gyre diversity in generated-world current-speed and SST
  fields under the current R2/R4 foundation, while preserving the phase4
  real-Earth improvements.

2026-07-06 - Replay-R acceptance protocol tightened to map-read attribution
Changed:
- Updated the Earth fitting plan and this climate plan to make direct visual
  reading of real-Earth reference and real-Earth replay maps a required phase
  gate.  Metrics are still required for regression, but they cannot promote a
  phase without a written spatial attribution of what the Earth maps show and
  which mechanism owns each residual.
- Corrected the workflow so R0-R8 tune one Earth subgraph at a time.  Generated
  worlds are not current fitting targets and should not drive parameter changes;
  they return only in R9 as promotion guardrails after the Earth replay subgraph
  is plausible.
Next:
- Continue by selecting the next Earth-only subgraph, rendering the real
  reference map and matching Earth replay map side by side, reading the visible
  residuals, and then changing only the mechanism that explains those residuals.
- Immediate target at this checkpoint was R2 seasonal 10 m wind vector/speed
  replay on real Earth against `earth__seasonal_wind_u10_v10`.  This was later
  superseded by the R2a/R2b split: first fit seasonal SLP / pressure-source
  geometry, then return to wind translation.  R3 OSCAR currents, R4 OISST SST,
  R5 moisture, R6 precipitation, R7 sea ice, and R8 classes/biomes stay blocked
  until their upstream Earth-only map-read packet passes.

2026-07-06 - Earth-only subgraph fitting contract restated
Changed:
- Moved the active fitting contract to the front of
  `docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md` so it is clear that R0-R8 are
  real-Earth replay phases, not virtual-world tuning phases.
- Renamed generated terminal-world inputs in the Earth fitting plan as R9-only
  guardrail inputs.  They are not allowed to choose R0-R8 mechanism changes.
- Updated this plan's Earth-based fitting track to state that a mechanism
  change during R0-R8 must be justified by the matching real-Earth reference
  subgraph, the Aevum real-Earth replay subgraph, and their residual/error maps.
- Added the current active R2 wind packet: seasonal 10 m wind speed/vector,
  eastward/northward components, speed delta, and vector-error maps against
  `earth__seasonal_wind_u10_v10`.
 Next:
- Continue only with the R2 Earth wind replay packet.  Do not use generated
  terminal worlds, precipitation maps, sea-ice maps, or biome maps as the basis
  for the next code change.

2026-07-06 - R2 Earth wind replay geography modulation checkpoint
Changed:
- Kept the active work inside R2 wind/pressure mechanics.  No generated worlds,
  precipitation thresholds, sea-ice thresholds, Koppen classes, or biome
  thresholds were used to justify this change.
- `_seasonal_circulation` now keeps the background wind belts continuous and
  moves the Southern Hemisphere near-surface westerly core poleward with a
  broader width, matching the real-Earth Southern Ocean wind belt better than a
  symmetric midlatitude center.
- `_geographic_circulation_anomalies` now receives geography primitives and
  adds a deterministic land/ocean stationary response: open oceans strengthen
  westerlies and trades, while continental interiors apply season-specific
  near-surface drag to the current wind vector.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_wind_seasonal_stationary_20260706/`.
- Wind diagnostic output:
  `out_real_earth_wind_replay_r2_seasonal_stationary_20260706/`.
- Direct map read: the hard latitude-band wall is reduced and the Southern
  Ocean band is more poleward, but the replay is still too stripe-like and land
  wind structure remains wrong.
- Wind metrics versus the continuous-asymmetric checkpoint: seasonal speed MAE
  `2.22 -> 1.99 m/s`, vector RMSE `3.68 -> 3.41 m/s`, direction-cosine p50
  `0.832 -> 0.853`, all-cell speed-pattern correlation `0.154 -> 0.327`, and
  ocean speed-pattern correlation `0.201 -> 0.368`.  Land speed-pattern
  correlation remains negative at `-0.154`, so R2 is not accepted.
- Real-Earth replay validation passes.  Temperature remains observer-only here:
  surface-temperature MAE is `2.63 C`, and the Earth-aware latitude jump remains
  below the same-grid Earth Antarctic/Southern-Ocean jump.  Annual wind p90
  remains low at `5.02/6.72 m/s` replay/Earth.
- Tests: `py_compile` passes for the edited module.  Target tests pass:
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`,
  `tests/test_real_earth_wind_replay.py`, and
  `tests/test_real_earth_climate_replay.py` (`5 passed in 90.69s`).
Next:
- Stay in R2.  Add or render explicit pressure-center/stationary-wave
  diagnostics against Earth pressure/wind structure, then repair the land
  stationary-wave residual.  Do not proceed to R3 OSCAR currents until the R2
  wind maps stop reading as mostly latitude bands.

2026-07-06 - R2 stationary-wave and pressure-center diagnostic
Changed:
- Extended `real-earth-wind-replay` with Earth-only R2 stationary-wave
  diagnostics: wind-speed zonal anomaly, eastward-wind zonal anomaly,
  standardized Earth SLP anomaly, replay pressure-proxy anomaly, pressure
  zonal anomaly, and corresponding residual maps.
- Added metrics for speed/eastward zonal-anomaly correlations and standardized
  pressure/pressure-zonal correlations.  The diagnostic remains optional for
  reference packages without SLP fields, but the R6 Earth package has the
  required seasonal SLP arrays.
- Tested a direct boundary-layer pressure-wave response inside R2 and rejected
  it because it worsened speed-zonal-anomaly correlation and vector RMSE.  That
  physical edit was reverted; the diagnostic extension was retained.
Evidence:
- R2 replay output at this diagnostic checkpoint:
  `out_real_earth_climate_replay_replay_r2_stationary_wave_current_20260706/`.
- Wind/pressure diagnostic output at this diagnostic checkpoint:
  `out_real_earth_wind_replay_r2_stationary_wave_current_20260706/`.
- Metrics: seasonal speed MAE `1.99 m/s`, vector RMSE `3.41 m/s`,
  direction-cosine p50 `0.853`, speed pattern correlation `0.327`, ocean speed
  pattern correlation `0.368`, and land speed pattern correlation `-0.154`.
  Pressure zonal-anomaly correlation is moderate (`0.561` all, `0.578` land,
  `0.547` ocean), while wind-speed zonal-anomaly correlation is weak (`0.148`
  all, `0.060` land, `0.155` ocean).
- Map read: Earth wind-speed zonal anomalies show broad continent/basin-scale
  patches; Aevum replay mostly shows weak coastal/topographic remnants and
  latitude-band texture.  Replay pressure centers are closer to Earth than the
  wind-speed stationary waves, so the next owner is the pressure/roughness/
  geostrophic-to-surface wind translation.
- Validation/tests: real-Earth climate replay passes validation.  Target tests
  pass: `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 93.06s`).  `py_compile` passes for the edited diagnostics and
  climate module.
Next:
- Stay in R2.  Design a proper pressure-gradient plus surface-roughness/
  boundary-layer wind-speed response and rerun the same stationary-wave
  diagnostic.  R3 currents remain blocked.

2026-07-06 - R2 scalar stationary roughness checkpoint
Changed:
- Added a conservative scalar wind-speed modulation inside R2
  `_geographic_circulation_anomalies`.  It preserves wind direction and uses
  pressure zonal-anomaly amplitude, open-ocean exposure, continent interiority,
  terrain roughness, wind gaps, and latitude-band context to slightly amplify
  or damp near-surface wind speed.
- The response is multiplied by `geography_strength`, so tiny-island and
  waterworld cases do not get a fake continental stationary wave.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_scalar_stationary_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_scalar_stationary_20260706/`.
- Metrics versus the prior stationary-wave diagnostic: seasonal speed MAE
  `1.986 -> 1.968 m/s`, vector RMSE `3.414 -> 3.386 m/s`, all-cell speed
  pattern correlation `0.327 -> 0.334`, speed-zonal-anomaly correlation
  `0.148 -> 0.183`, land speed-zonal-anomaly correlation `0.060 -> 0.081`,
  and eastward-zonal-anomaly correlation `0.291 -> 0.304`.  Direction-cosine
  p50 stays effectively unchanged (`0.853`).
- Remaining blockers: land speed-pattern correlation remains negative
  (`-0.153`), ocean speed-zonal-anomaly correlation does not improve
  (`0.153`), and seasonal p90 wind speed remains low (`5.51` replay versus
  `7.00 m/s` Earth).  The map still reads too stripe-like, so this is an
  accepted checkpoint but not R2 promotion.
- Validation/tests: real-Earth climate replay passes validation.  Target tests
  pass: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`5 passed in 91.65s`).
Next:
- Stay in R2.  The next repair should address why pressure-center structure
  does not create sufficiently broad land/ocean wind-speed stationary waves,
  probably by replacing the current very local roughness response with a
  basin/continent-scale boundary-layer response.  R3 remains blocked.

2026-07-06 - R2 basin/continent boundary-layer checkpoint
Changed:
- Added a broader R2 wind-speed response on top of the scalar stationary
  roughness checkpoint.  The new term smooths pressure zonal-anomaly amplitude
  separately over land and ocean domains, then uses that basin/continent-scale
  support to modulate wind speed while preserving wind direction.
- Added a small open-ocean strong-tail support so improving stationary-wave
  structure does not further suppress the already-low p90 wind speed.
- The new term remains scaled by `geography_strength`, so tiny-island and
  waterworld cases stay protected from false continental stationary waves.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_basin_boundary_layer_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_basin_boundary_layer_20260706/`.
- Metrics versus the scalar stationary checkpoint: seasonal speed MAE
  `1.968 -> 1.945 m/s`, vector RMSE `3.386 -> 3.367 m/s`, all-cell
  speed-pattern correlation `0.334 -> 0.355`, land speed-pattern correlation
  `-0.153 -> -0.150`, speed-zonal-anomaly correlation `0.183 -> 0.236`, land
  speed-zonal-anomaly correlation `0.081 -> 0.113`, ocean speed-zonal-anomaly
  correlation `0.153 -> 0.172`, eastward-zonal-anomaly correlation
  `0.304 -> 0.339`, and seasonal p90 replay wind speed `5.51 -> 5.57 m/s`
  versus Earth `7.00 m/s`.
- Map read: replay wind-speed stationary waves remain too weak and too
  latitude-banded, but the broad response is stronger than the scalar-only
  checkpoint without adding visible noise or breaking the Southern Ocean band.
  R2 is still not promoted.
- Observer fields remain stable: real-Earth replay validation passes, annual
  wind p90 improves to `5.12/6.72 m/s` replay/Earth, and ocean-current p90
  ratio remains near `0.997`.
- Tests: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  pass (`5 passed in 90.94s`).
Next:
- Stay in R2.  Further repair should target the still-negative land speed
  pattern correlation and low strong-wind tail, ideally through continent-scale
  roughness/pressure placement rather than local vector pressure kicks.  R3
  remains blocked.

2026-07-06 - Replay-R plan narrowed to one Earth subgraph at a time
Changed:
- Restated the active fitting workflow in
  `docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md` as a single-subgraph Earth replay
  protocol.  At this checkpoint the active packet was still R2 seasonal 10 m
  wind vector/speed on the real-Earth grid against
  `earth__seasonal_wind_u10_v10`; this was later narrowed to R2a pressure-source
  geometry before R2b wind translation.
- Made the evidence order explicit: same-grid Earth reference map, same-grid
  Aevum real-Earth replay map, residual/vector-error maps, written map-read
  attribution, then metrics as screening/regression evidence.  Global means,
  scalar envelopes, or pass/fail tables cannot promote a packet whose maps
  still have the wrong geography.
- Clarified that generated terminal worlds are R9 guardrails only.  They must
  not choose R0-R8 parameters and should not be rendered during the active R2
  fitting loop unless a later R9 promotion check is explicitly being run.
- Reaffirmed the current R2 blockers: replay wind-speed stationary waves are
  too weak and too latitude-banded, land wind-speed pattern correlation remains
  negative, and the strong-wind tail remains low versus Earth.
Next:
- Continue only with R2 pressure/roughness/geostrophic-to-surface wind
  translation on real Earth.  Do not move to OSCAR currents, OISST SST,
  moisture, precipitation, sea ice, Koppen, biomes, or generated-world
  evaluation until the R2 Earth map-read packet is accepted.

2026-07-06 - R2 ocean-tail plus land-drag checkpoint, not promoted
Changed:
- Kept the active work inside the single R2 Earth wind packet.  No generated
  worlds, current/SST tuning, precipitation, sea ice, Koppen, or biome
  thresholds were used.
- Added a scalar near-surface wind-speed response in
  `_geographic_circulation_anomalies`: stronger open-ocean midlatitude
  storm-track tail plus stronger continent-interior/terrain roughness
  boundary-layer drag.  The response preserves wind direction and replaces no
  downstream mechanism.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_ocean_tail_land_drag_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_ocean_tail_land_drag_20260706/`.
- Metrics versus the basin-boundary-layer checkpoint: seasonal speed MAE
  `1.945 -> 1.935 m/s`, vector RMSE `3.367 -> 3.363 m/s`, seasonal replay p90
  `5.57 -> 5.90 m/s` versus Earth `7.00 m/s`, all-cell speed-pattern
  correlation `0.355 -> 0.401`, land speed-pattern correlation
  `-0.150 -> -0.127`, speed-zonal-anomaly correlation `0.236 -> 0.317`, land
  speed-zonal-anomaly correlation `0.113 -> 0.159`, ocean speed-zonal-anomaly
  correlation `0.172 -> 0.198`, and eastward-zonal-anomaly correlation
  `0.339 -> 0.377`.
- Map read: ocean wind-tail amplitude and continent weak-wind patches are less
  underpowered, but replay still reads too much like latitude belts compared
  with Earth's North Atlantic/North Pacific/Southern Ocean storm-track
  geography and continental stationary-wave structure.
- Tests: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  pass (`5 passed in 91.85s`).
Next:
- Stay in R2.  The remaining owner is not p90 amplitude alone; it is
  continent/basin-scale stationary-wave placement and wind-direction/vector
  structure.  R3 currents and all downstream maps remain blocked.

2026-07-06 - R2 broad pressure-steering checkpoint, not promoted
Changed:
- Kept the active work inside the R2 Earth wind packet.  Added a small
  broad-scale geostrophic steering response from seasonal pressure zonal
  anomalies inside `_geographic_circulation_anomalies`.
- The steering term is Coriolis-damped near the equator, supported by
  trade/westerly belts plus land/ocean exposure, and capped through the
  existing wind-speed envelope.  It is intended to bend the latitude-band wind
  field around continent/basin pressure centers without adding the rejected
  local pressure-vector kick.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_pressure_steering_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_pressure_steering_20260706/`.
- Metrics versus the ocean-tail plus land-drag checkpoint: speed MAE
  `1.935 -> 1.932 m/s`, vector RMSE `3.363 -> 3.354 m/s`,
  direction-cosine p50 `0.853 -> 0.856`, direction-cosine p10
  `-0.673 -> -0.667`, seasonal replay p90 `5.90 -> 5.92 m/s` versus Earth
  `7.00 m/s`, land speed-pattern correlation `-0.127 -> -0.112`,
  speed-zonal-anomaly correlation `0.317 -> 0.323`, and
  eastward-zonal-anomaly correlation `0.377 -> 0.382`.
- Map read: no obvious new noise or wrong local pressure spirals, and
  component residuals improve slightly.  The replay still reads too
  latitude-banded relative to Earth's North Atlantic/North Pacific/Southern
  Ocean storm tracks and continental stationary-wave structure, so R2 is not
  promoted.
- Tests: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  pass (`5 passed in 92.70s`).
Next:
- Stay in R2.  Further progress probably needs better continent/basin pressure
  placement or stationary-wave source geometry, not stronger scalar wind-speed
  amplification.  R3 currents remain blocked.

2026-07-06 - R2 solstice basin-pressure source checkpoint, not promoted
Changed:
- Stayed inside the R2 Earth wind/pressure packet.  Added a geography-derived
  ocean pressure source: large open ocean basins support solstice-season
  subpolar lows in the winter hemisphere and broad subtropical ocean highs.
- The source is derived from basin id, coast distance, latitude band, and
  season; it is not fitted from generated worlds and does not use downstream
  current/SST/precipitation/biome compensation.  A shoulder-season version was
  tested and rejected because it degraded MAM/SON pressure placement.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_solstice_basin_pressure_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_solstice_basin_pressure_20260706/`.
- Metrics versus the previous pressure-steering checkpoint: speed MAE
  `1.932 -> 1.927 m/s`, vector RMSE `3.354 -> 3.343 m/s`,
  direction-cosine p10 `-0.667 -> -0.662`, seasonal replay p90 remains about
  `5.92 m/s`, land speed-pattern correlation `-0.112 -> -0.110`,
  speed-zonal-anomaly correlation `0.323 -> 0.326`, ocean speed-zonal-anomaly
  correlation `0.208 -> 0.211`, pressure standardized ocean correlation
  `0.083 -> 0.135`, pressure-zonal ocean correlation `0.546 -> 0.560`,
  DJF pressure-zonal correlation `0.668 -> 0.682`, and JJA pressure-zonal
  correlation `0.712 -> 0.722`.
- Map read: the upstream pressure owner is improved in the solstice seasons,
  and no visible local noise was added.  Replay wind maps still read too much
  like latitude belts; land speed-pattern correlation remains negative, and
  basin/storm-track placement is still not Earthlike enough for R2 promotion.
- Tests: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  pass (`5 passed in 92.35s`).
Next:
- Stay in R2.  The next owner is storm-track/stationary-wave geometry and
  land negative-correlation attribution, not downstream ocean-current tuning.

2026-07-06 - R2 land boundary-layer and polar katabatic checkpoint, not promoted
Changed:
- Stayed inside the R2 Earth wind/pressure packet.  Component-level residual
  attribution showed that non-polar land was too windy while Antarctica was
  strongly under-windy.
- Added an extra non-polar land boundary-layer drag and a polar highland
  downslope katabatic wind term.  The katabatic support is derived from
  latitude, smoothed elevation, terrain gradient, and continent interiority; it
  is not a generated-world or downstream climate-class compensation.
Evidence:
- Earth replay output:
  `out_real_earth_climate_replay_replay_r2_land_katabatic_v2_20260706/`.
- Wind/pressure diagnostic output:
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/`.
- Metrics versus the solstice basin-pressure checkpoint: speed MAE
  `1.927 -> 1.865 m/s`, vector RMSE `3.343 -> 3.289 m/s`,
  direction-cosine p10 `-0.662 -> -0.654`, seasonal replay p90 remains about
  `5.91 m/s` versus Earth `7.00 m/s`, all-cell speed-pattern correlation
  `0.405 -> 0.444`, land speed-pattern correlation `-0.110 -> 0.145`,
  speed-zonal-anomaly correlation `0.326 -> 0.367`, land speed-zonal-anomaly
  correlation `0.172 -> 0.210`, ocean speed-zonal-anomaly correlation
  `0.211 -> 0.227`, and eastward-zonal-anomaly correlation `0.384 -> 0.416`.
- Map read: replay now expresses Antarctic/polar highland wind support and
  non-polar land drag more plausibly, with no visible new local noise.  R2 is
  still not promoted because the maps remain too latitude-banded overall and
  ocean storm-track/basin-scale stationary waves are still weaker than Earth.
- Tests: `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_engine.py::test_seasonal_winds_migrate_itcz_and_storm_tracks`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  pass (`5 passed in 92.69s`).
Next:
- Stay in R2.  The next likely owner is ocean storm-track/basin-scale
  stationary-wave geometry; R3 currents remain blocked.

2026-07-06 - R2 ocean storm-track residual attribution, no promotion
Changed:
- No additional physics change was accepted after the land/katabatic
  checkpoint.  Two ocean storm-track scalar probes were rejected.
Evidence:
- Sector attribution from
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/` shows the ocean
  residual is a spatial-placement problem.  Near/coastal ocean is generally too
  weak, while some non-leeward midlatitude open-ocean belts are too strong.
- A land-leeward ocean corridor boost inferred from current wind direction
  worsened ocean speed-pattern, ocean zonal-anomaly, and eastward-anomaly
  correlations despite sometimes lowering scalar MAE.
- A near-coast boost plus open-ocean damping probe had the same failure mode:
  better scalar envelopes in some settings, but worse ocean geography.
Next:
- Stay in R2.  The next accepted change should improve storm-track source
  geometry itself, probably through baroclinic/coastal thermal-front and
  basin-sector support derived from real Earth replay fields.  Do not proceed
  to R3 currents or compensate with downstream maps.

2026-07-06 - Earth-only single-subgraph workflow tightened
Changed:
- Updated the front of this plan to make real-Earth single-subgraph replay the
  authoritative workflow for R0-R8.  Each loop must begin from one same-grid
  Earth reference map, the matching Aevum replay map, and residual/error maps,
  then record a map-read attribution before code changes.
- Clarified that generated worlds, old C5 gates, global means, p50/p90
  envelopes, correlations, and pass/fail tables are regression evidence only.
  They may reject a change, but they cannot promote a visibly wrong Earth
  subgraph.
- Updated the current R2 packet: land speed-pattern correlation has turned
  positive at the land/katabatic checkpoint, so the active blocker is now
  ocean storm-track and basin-scale stationary-wave placement rather than
  generic land negative-correlation repair.
Next:
- Stay in R2 and fit the real-Earth seasonal 10 m wind subgraph only.  The next
  accepted mechanism must target storm-track source geometry and
  pressure-to-surface-wind translation on the real-Earth grid.  R3 currents,
  R4 SST, R5 moisture, R6 precipitation, sea ice, Koppen, biomes, and generated
  worlds remain blocked as tuning targets.

2026-07-06 - R2 storm-track post-processing probes rejected
Changed:
- No model code was accepted.  The loop stayed inside the real-Earth R2
  pressure/wind packet and tested candidate fixes offline against the current
  land/katabatic checkpoint.
Evidence:
- Current checkpoint remains
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/`, paired with
  `out_real_earth_climate_replay_replay_r2_land_katabatic_v2_20260706/`.
- Baroclinic/thermal-front scalar modifiers worsened ocean speed-pattern,
  ocean speed-zonal, and ocean eastward-anomaly correlations even at low
  amplitude.  They are rejected as another misplaced multiplier.
- Extra ocean pressure-steering slightly improved ocean anomaly metrics but did
  not visibly change the latitude-banded map read; probe output:
  `out_probe_real_earth_wind_replay_r2_pressure_steering_ocean_extra2_20260706/`.
- East-coast winter ocean-low pressure sources did not improve the real-Earth
  pressure map against `earth__seasonal_slp_anomaly_hPa`.
- Orographic downstream stationary-wave redistribution plus modest ocean tail
  improved some scalar pattern correlations but still rendered as latitude
  belts and degraded eastward anomaly structure; probe output:
  `out_probe_real_earth_wind_replay_r2_orographic_wave_tail_20260706/`.
- Broad open-ocean tail amplification can move p90 wind speed toward Earth, but
  it worsens MAE and eastward-anomaly placement as amplitude grows.
Next:
- Split the active R2 work into R2a pressure/source geometry followed by R2b
  wind translation.  Next fit the real-Earth SLP/pressure-source map directly
  before accepting more wind-speed or vector modifiers.  R3 and all downstream
  maps remain blocked as tuning targets.

2026-07-06 - R2a pressure/source diagnostic checkpoint
Changed:
- Extended `real-earth-wind-replay` with pressure-only R2a evidence:
  `pressure_standardized_delta_seasons.png`,
  `real_earth_pressure_replay_contact_sheet.png`, and standardized-pressure
  MAE metrics for all cells, land, ocean, and each season.
- Updated `tests/test_real_earth_wind_replay.py` to require the new pressure
  residual assets.
Evidence:
- Current R2a baseline output:
  `out_real_earth_pressure_replay_r2a_current_20260706/`.
- Map read: replay pressure captures the broad seasonal land thermal signal,
  but remains too smooth and too blocky.  It lacks enough ocean basin pressure
  centers and mountain/coast-driven stationary-wave structure, which explains
  why the wind replay still reads as latitude belts.
- Current metrics: pressure standardized MAE all/ocean `0.317/0.277`,
  standardized correlation all/land/ocean `0.539/0.668/0.135`, pressure
  zonal-anomaly correlation all/land/ocean `0.567/0.578/0.560`, and seasonal
  pressure-zonal correlations `0.682/0.563/0.722/0.388` for
  DJF/MAM/JJA/SON.
- Pressure-source probes were rejected: stronger open-ocean winter subpolar
  lows improved ocean standardized correlation only to `0.160` while worsening
  MAE and leaving the map visually similar; reduced land thermal smoothing plus
  stronger ocean lows could raise ocean standardized correlation to about
  `0.171`, but degraded all-cell/land/zonal structure; basin-sector
  downstream plumes, mountain downstream waves, and ocean thermal-front lows
  worsened the pressure map.
- Tests: `tests/test_real_earth_wind_replay.py` and
  `tests/test_real_earth_climate_replay.py` pass (`2 passed in 4.95s`).
Next:
- Stay in R2a.  The next accepted mechanism should add explicit basin-scale
  pressure-center or stationary-wave pressure objects derived from Earth
  geography inputs, then be judged by the pressure contact sheet before any
  R2b wind translation or downstream tuning.

2026-07-06 - R2a active packet restated after workflow correction
Changed:
- Tightened both climate planning documents so the current task is no longer
  described as a mixed wind/pressure fitting loop.  The active target is exactly
  one real-Earth subgraph: R2a seasonal SLP / pressure-source geometry.
- Reclassified the best wind checkpoint
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/` as downstream
  observer evidence and historical motivation for moving upstream, not as the
  current fitting target.
- Made R2b wind translation, R3 currents, R4 SST, R5 moisture, R6
  precipitation, sea ice, Koppen, biomes, and generated worlds explicitly
  blocked as tuning targets until the R2a pressure contact sheet passes map
  read.
Next:
- Read the real-Earth seasonal SLP anomaly map, the replay pressure-proxy map,
  and the pressure residual map first.  The next code step may expose
  pressure-center / stationary-wave objects for diagnosis, but any accepted
  physics change must be judged only against the R2a pressure maps.

2026-07-07 - R2a pressure-center diagnostic object checkpoint
Changed:
- Added diagnostic-only R2a pressure/source geometry fields to
  `ClimateModule`: `atmosphere.pressure_center_support`,
  `atmosphere.pressure_center_id`, and
  `atmosphere.stationary_wave_pressure_support`.
- Added `atmosphere.pressure_centers` objects with season, high/low kind,
  land/ocean domain, centroid, area, dominant basin/continent, and coast/barrier
  summaries.
- Archived the new fields in `terminal_climate_arrays.npz`, wrote
  `pressure_centers.json`, and extended `real-earth-wind-replay` so the R2a
  pressure contact sheet can include pressure-center and stationary-wave support
  panels.
Evidence:
- Targeted tests pass:
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_real_earth_climate_replay.py`, and
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`.
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_pressure_centers_20260706/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_pressure_centers_20260706/`.
- The replay emits 45 pressure-center objects: 25 lows, 20 highs; 27 land, 17
  ocean, and 1 mixed.
Map read:
- The diagnostic support panels are useful for attribution, but they confirm
  that R2a pressure is not yet accepted.  Replay pressure still reads as broad
  continental blobs plus strong latitude bands, especially near the Southern
  Ocean.  Compared with Earth SLP, ocean basin pressure centers and
  coast/mountain stationary-wave structure remain too weak and too smooth.
Next:
- Stay in R2a.  Use the new pressure-center and stationary-wave objects to
  design the first accepted pressure-source geometry change.  Do not tune R2b
  wind translation, currents, SST, moisture, precipitation, sea ice, Koppen,
  biomes, or generated worlds until the pressure contact sheet itself improves.

2026-07-07 - R2a M2 pressure-genesis v6 checkpoint, not promoted
Changed:
- Kept the active work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.  No generated worlds and no downstream tuning were used as
  acceptance evidence.
- Added auditable M2 causal source archives:
  `atmosphere.pressure_genesis_source`,
  `atmosphere.ocean_pressure_low_source_support`,
  `atmosphere.ocean_pressure_high_source_support`,
  `atmosphere.land_pressure_source_support`, and
  `atmosphere.terrain_pressure_wave_source_support`.
- Added bounded source-to-pressure projection archived as
  `atmosphere.pressure_genesis_wave_transfer`.  The transfer is deliberately
  gated over polar ice-cap edge geometry to avoid turning Antarctic margin
  artifacts into ordinary continental stationary waves.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v6_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Source and transfer maps:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_pressure_genesis_source_seasons.png` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- Major ocean semantic basin support and M1 energy support are now readable.
  DJF North Pacific / North Atlantic low-source patches and JJA Southern Ocean
  source sectors are visible.  The v6 transfer reduces the prior Antarctic
  polar-edge artifact.
- The final pressure proxy still remains too smooth compared with real Earth
  SLP.  Continental transfer is still too block-like, and terrain/coast/SST
  front waveguides do not yet organize Aleutian, Icelandic, Southern Ocean, and
  continental pressure centers strongly enough.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.318/0.419/0.277`.
- Standardized-pressure correlation all/land/ocean:
  `0.548/0.670/0.209`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.579/0.589/0.570`.
- Wind observer metrics remain downstream-only: seasonal speed MAE
  `1.866 m/s`, direction cosine p50 `0.858`.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.47s`).
Next:
- R2a is still not promoted.  The next accepted change should stay in M2 and
  improve source-to-pressure spatial propagation using terrain, coast, SST
  front, and basin waveguides.  R2b wind, R3 currents, R4 SST, R5 moisture, R6
  precipitation, R7 feedback, R8 classes/biomes, and generated worlds remain
  blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v7 waveguide checkpoint, not promoted
Changed:
- Kept the active work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a local weighted-diffusion helper for masked fields and used it only
  inside M2 pressure genesis.
- Replaced the v6 source-to-pressure transfer with a geography-weighted
  waveguide transfer.  Ocean propagation is weighted by SST-front support,
  same-latitude SST anomaly, open-ocean exposure, and subpolar/basin support.
  Land propagation is weighted by coast-strength, terrain barriers,
  land-source support, and terrain-wave support, with the existing non-polar
  gate.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v7_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- The v7 transfer is more organized along Northern Hemisphere coastal/storm
  track waveguides and does not reintroduce the Antarctic polar-edge artifact.
  The source map remains auditable and unchanged in intent from v6.
- The final pressure map is still too similar to v6 and too smooth/block-like
  versus Earth SLP.  R2a is therefore not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.317/0.418/0.277`.
- Standardized-pressure correlation all/land/ocean:
  `0.549/0.671/0.209`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.581/0.591/0.571`.
- Wind observer metrics remain downstream-only: seasonal speed MAE
  `1.866 m/s`, direction cosine p50 `0.858`.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.79s`).
Next:
- Stay in R2a/M2.  The next change should build object-level pressure-center
  projection on top of the v7 waveguide baseline so Aleutian, Icelandic,
  Southern Ocean, and continental thermal centers become visible in the final
  pressure map.  Do not tune R2b/R3/R4/R5/R6/R7/R8 or generated worlds.

2026-07-07 - R2a M2 pressure-genesis v8 object-projection checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added bounded object-level pressure-center projection on top of the v7
  weighted waveguide transfer.  The projection uses existing causal supports:
  ocean-low support, land thermal-center source support, and terrain-wave
  source support.
- Increased the effective M2 source-to-pressure transfer enough for projected
  objects to affect the final pressure proxy, while keeping clipping and the
  non-polar land gate.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v8_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF North Pacific / North Atlantic low-pressure footprints and Eurasian
  winter-high transfer are clearer than v7.  JJA Southern Ocean remains
  segmented and wave-like; no Antarctic edge artifact returned.
- Final standardized pressure is still too smooth versus Earth SLP, so R2a is
  not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.317/0.418/0.277`.
- Standardized-pressure correlation all/land/ocean:
  `0.550/0.672/0.212`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.581/0.592/0.570`.
- Wind observer metrics remain downstream-only: seasonal speed MAE
  `1.866 m/s`, direction cosine p50 `0.858`.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 70.93s`).
Next:
- Stay in R2a/M2.  The next change should improve pressure-center placement
  and anisotropic footprint geometry using source object centroids,
  coast/terrain orientation, SST-front orientation, and basin/continent labels.
  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v9 directional-footprint checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added directional masked diffusion helpers for M2 pressure genesis.
- Replaced v8's ordinary weighted diffusion in source-to-pressure transfer
  with directional diffusion.  Ocean transfer follows a storm-track/SST-front
  axis; land transfer follows terrain-barrier and coastal axes.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v9_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- v9 preserves v8's DJF North Pacific / North Atlantic low-pressure footprints
  and Eurasian winter-high transfer.  JJA Southern Ocean remains segmented and
  no Antarctic polar-edge artifact returned.
- Final standardized pressure still looks too similar to v8 and too smooth
  versus Earth SLP, so R2a is not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.317/0.418/0.277`.
- Standardized-pressure correlation all/land/ocean:
  `0.550/0.672/0.212`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.581/0.592/0.570`.
- Wind observer metrics remain downstream-only: seasonal speed MAE
  `1.866 m/s`, direction cosine p50 `0.858`.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.57s`).
Next:
- Stay in R2a/M2.  Anisotropy is now present but too weak to solve the final
  pressure map.  The next change should improve final-pressure expression of
  source-supported pressure-center objects without global amplitude painting.
  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v10 nonzonal-expression checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a bounded final-pressure expression solve after the v9 directional
  transfer.  The solve enhances source-supported nonzonal pressure anomaly only
  where the anomaly sign aligns with the signed M2 source/transfer field, then
  removes the latitude-band mean and clips the increment.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v10_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF North Pacific / North Atlantic lows and the Eurasian winter-high transfer
  are more forceful than v9.  JJA Southern Ocean remains segmented and
  wave-like, with no Antarctic edge artifact.
- Final standardized pressure improves but remains too smooth relative to
  Earth SLP, so R2a is not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.317/0.416/0.276`.
- Standardized-pressure correlation all/land/ocean:
  `0.550/0.671/0.220`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.582/0.593/0.571`.
- Wind observer metrics remain downstream-only: seasonal speed MAE
  `1.866 m/s`, direction cosine p50 `0.858`.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.19s`).
Next:
- Stay in R2a/M2.  v10 becomes the current baseline because it improves MAE,
  ocean pressure correlation, and zonal-anomaly correlation together.  The next
  change should improve pressure-center placement and shape in the final
  pressure map, with a guard against high-latitude expression becoming too
  band-like.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as
  fitting targets.

2026-07-07 - R2a M2 pressure-genesis v16 downwind-cold-continent checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- v11/v12 tried object-local final-pressure expression and made the map no
  more readable while degrading pressure metrics; v13/v14 tried stronger object
  projection/source reinforcement and also failed promotion.
- v16 keeps the v10 final-pressure expression solve and adds one mechanism in
  M2 source generation: winter cold continental highs are spread downwind along
  the westerly/storm-track axis into adjacent open ocean, then used only as
  support for subpolar ocean-low source selection.  Southern Hemisphere
  downwind support is damped so the Southern Ocean front/wavenumber gate remains
  the primary Southern Ocean control.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v16_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF North Pacific / North Atlantic source placement remains readable and is
  now partly derived from cold-continent downstream support rather than only
  open-ocean latitude/front scoring.
- No new Antarctic edge artifact or full Southern Ocean ring was introduced.
- Final standardized pressure remains too smooth relative to Earth SLP, so R2a
  is not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.317/0.416/0.276`.
- Standardized-pressure correlation all/land/ocean:
  `0.551/0.671/0.221`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.583/0.594/0.572`.
- Compared with v10, v16 is essentially neutral on all-MAE, slightly improves
  ocean MAE and pressure correlations, and improves JJA zonal-anomaly
  correlation, but the visual gap is still not closed.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.24s`).
Next:
- Stay in R2a/M2.  The remaining owner is now explicit: the M2 source maps are
  readable, but the final pressure proxy is still dominated by the smoothed
  upstream pressure field.  The next accepted change should re-balance M2
  source-supported pressure against the upstream smoothed proxy without broad
  latitude-band painting.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain
  blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v17 source/transfer rebalance checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Used v16 arrays to reconstruct `final = upstream_base + M2_source +
  M2_transfer` and probe the source/transfer balance before editing code.
- The probe showed that strengthening raw source support directly degrades the
  SLP pattern.  v17 therefore keeps source support as the causal trigger,
  reduces its direct pressure contribution to `0.80`, and increases the
  bounded source-to-pressure transfer expression to `1.45`.  The archived
  `atmosphere.pressure_genesis_wave_transfer` now reflects this stronger
  transfer contribution.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v17_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- The transfer map is visibly stronger and more active in the already
  source-supported Northern Hemisphere subpolar regions and the JJA Southern
  Ocean wave belt.
- Final pressure is incrementally less dominated by the smooth upstream proxy,
  especially in JJA, but it remains too smooth and too weakly center-organized
  relative to Earth SLP.  R2a is therefore still not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.316/0.416/0.275`.
- Standardized-pressure correlation all/land/ocean:
  `0.552/0.672/0.222`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.584/0.595/0.573`.
- Compared with v16, v17 improves all-pressure MAE, ocean MAE, all/land/ocean
  pressure correlations, and all/land/ocean zonal-anomaly correlations.  DJF
  and SON zonal-anomaly correlations are slightly lower, so this remains a
  checkpoint rather than an acceptance point.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 104.60s`).
Next:
- Stay in R2a/M2.  v17 confirms that M2 transfer should carry more of the
  source-to-pressure expression, but final pressure centers still need better
  coherence and placement.  The next accepted change should improve pressure
  center morphology without introducing hard-coded seasonal gain tables or
  broad latitude-band painting.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds
  remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v18 transfer-morphology checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Used v17 arrays to probe M2-only morphology terms.  The retained direction is
  a bounded adjustment to the already-projected transfer: slightly deepen
  ocean-low cores where low support and negative transfer overlap, while
  damping broad subtropical-high, land-core, and terrain-wave transfer
  contributions that make final pressure too smooth or fragmented.
- No hard-coded coordinates, no hard-coded seasonal gain table, and no R2b/R3/
  SST/precipitation/biome tuning was introduced.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v18_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- The transfer panel is a little more center-like and less dominated by terrain
  fragments.  No new Southern Ocean ring, Antarctic edge artifact, or speckled
  high-latitude band was introduced.
- Final pressure is modestly more coherent than v17, but still smoother and
  less center-organized than Earth SLP.  R2a is therefore still not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.316/0.416/0.275`.
- Standardized-pressure correlation all/land/ocean:
  `0.553/0.673/0.224`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.585/0.596/0.574`.
- Compared with v17, v18 improves all/land MAE and all/land/ocean pressure
  and zonal-anomaly correlations.  Ocean MAE is nearly unchanged but slightly
  worse, and JJA zonal-anomaly correlation steps back from the v17 high while
  staying above v16.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 204.10s`).
Next:
- Stay in R2a/M2.  The next owner is still pressure-center morphology and
  placement in the final pressure proxy.  v18 is a balanced checkpoint, but R2a
  acceptance still requires Earthlike Aleutian/Icelandic/Southern Ocean and
  continental seasonal centers to be visually coherent in the final SLP replay.
  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v19 thermal-phase checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a bounded thermal-phase adjustment to the already-projected M2 transfer:
  same-latitude SST and land-temperature anomalies now modulate whether
  high-support, low-support, and land-source footprints strengthen or weaken
  the final pressure expression.  Cold subtropical ocean and cold continental
  interiors can support highs; warm ocean/land phases avoid promoting the same
  support into the wrong seasonal sign.
- This remains M2-only.  No coordinate-specific target, hard-coded seasonal
  gain table, R2b/R3/SST/precipitation/biome tuning, or generated-world fitting
  was introduced.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v19_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- v19 improves seasonal differentiation in the transfer and final pressure
  panels, especially MAM/SON phase behavior, without introducing a full
  Southern Ocean ring, Antarctic edge artifact, or heat-wall latitude band.
- Final pressure is still too smooth and not center-organized enough compared
  with Earth SLP, so R2a is still not promoted.  Wider subtropical
  thermal-phase patches in the transfer panel are a watch item for the next
  morphology pass.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.310/0.416/0.268`.
- Standardized-pressure correlation all/land/ocean:
  `0.564/0.675/0.269`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.587/0.600/0.574`.
- Compared with v18, v19 improves all/ocean MAE and all/land/ocean pressure
  correlation, with the largest gain in ocean standardized-pressure
  correlation.  MAM pressure MAE and zonal-anomaly correlation improve
  materially; DJF/JJA/SON zonal-anomaly correlations are nearly flat to
  slightly lower.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 139.05s`).
Next:
- Stay in R2a/M2.  v19 is now the current worktree checkpoint, but not an
  acceptance point.  The next owner remains final-pressure center coherence and
  placement: Aleutian/Icelandic/Southern Ocean and continental seasonal centers
  need to become coherent in the final SLP replay, not merely visible in
  source/transfer support maps.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds
  remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v25 domain-weighted source checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Tightened the Southern Ocean wavenumber/front gate by lowering its baseline
  floor and sharpening the wave core.  This makes the Southern Ocean source
  less ring-like in JJA while preserving the sector/front basis.
- Changed final pressure expression from a uniform direct-source weight to a
  domain-weighted source expression: land and Southern Ocean source expression
  are reduced, while North Hemisphere mid/high-latitude ocean-low cores retain
  more direct expression.  The intent is to keep Aleutian/Icelandic-type lows
  from disappearing while preventing broad land and Southern Ocean source
  fields from overprinting the final SLP proxy.
- Rejected intermediate trials: v20/v21 over-activated MAM/SON ocean lows and
  produced an unrealistic shoulder-season Southern Ocean/North Atlantic low
  field; v23 basin-area reweighting gave no useful visual gain and slightly
  hurt ocean correlations.
- This remains M2-only.  No coordinate-specific target, hard-coded seasonal
  gain table, R2b/R3/SST/precipitation/biome tuning, or generated-world fitting
  was introduced.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v25_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- JJA Southern Ocean no longer receives as much direct low-source overprint in
  the final pressure proxy; residuals there are smaller than v19.
- DJF Aleutian remains close to v19, so the domain weighting avoids the broad
  degradation caused by the uniform v24 source reduction.  DJF Icelandic Low is
  still too weak, and SON North Pacific still lacks a meaningful source.
- Final pressure remains smoother and less center-organized than Earth SLP, so
  v25 is retained only as a checkpoint.  R2a is still not promoted.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.309/0.416/0.266`.
- Standardized-pressure correlation all/land/ocean:
  `0.564/0.674/0.272`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.589/0.600/0.577`.
- Compared with v19, v25 improves all/ocean MAE, all/ocean standardized
  pressure correlation, and all/land/ocean zonal-anomaly correlation.  Land
  MAE/correlation are nearly flat to slightly worse, and DJF zonal-anomaly
  correlation is slightly lower.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 115.70s`).
Next:
- Stay in R2a/M2.  v25 is the current worktree checkpoint, but not an
  acceptance point.  The next owner remains missing/weak pressure-source
  placement, especially DJF Icelandic Low strength and SON North Pacific source
  emergence, plus final-pressure center compactness.  R2b/R3/R4/R5/R6/R7/R8
  and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v30 shoulder-season source-placement checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a winter Northern Hemisphere subpolar coastal/front supplement after
  low-support normalization.  This addresses the v25 attribution that the
  Icelandic Low was suppressed because the source used `open_ocean` too
  strongly in a basin-edge/frontal region.
- Added a bounded shoulder-season warm-ocean/front low-source candidate.  It
  uses positive seasonal SST anomaly, SST-front support, open-ocean exposure,
  subpolar latitude, and basin component selection, then applies a small
  post-object amplitude instead of letting the candidate normalize to a winter
  low.  This creates a weak SON North Pacific/North Atlantic low source without
  reactivating MAM lows.
- Rejected intermediate trials: v26 improved Icelandic Low but weakened
  Aleutian through global low-support normalization; v28 made SON ocean lows
  far too deep; v29 was closer but still too strong.  v30 keeps the same
  mechanism with lower shoulder-season amplitude.
- This remains M2-only.  No coordinate-specific target, hard-coded seasonal
  gain table, R2b/R3/SST/precipitation/biome tuning, or generated-world fitting
  was introduced.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v30_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- SON North Pacific now has a weak, local low-source/transfer expression; v25
  had essentially none.  MAM North Pacific/North Atlantic remain untriggered,
  so the shoulder-season mechanism is not simply painting both equinoxes.
- DJF Icelandic Low is slightly stronger than v25 while DJF Aleutian is not
  visibly degraded.  The final pressure field remains too smooth, and the
  Icelandic center is still too weak compared with Earth.
- The tradeoff is explicit: SON/ocean scalar correlations step back from v25,
  even though the missing SON source object is now present.  v30 is therefore a
  source-placement checkpoint, not an R2a acceptance point.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.309/0.416/0.266`.
- Standardized-pressure correlation all/land/ocean:
  `0.563/0.674/0.268`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.585/0.598/0.569`.
- Compared with v25, v30 improves DJF MAE slightly and introduces the missing
  SON North Pacific source, but all/ocean pressure correlations and SON
  zonal-anomaly correlation regress.  The regression is accepted only as a
  temporary mechanism checkpoint because the visual source-placement blocker
  moved in the right direction.
Tests:
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.16s`).
Next:
- Stay in R2a/M2.  v30 is the current worktree checkpoint, but not an
  acceptance point.  The next owner is final-pressure expression and center
  compactness: retain the newly present SON source while recovering ocean/zonal
  correlations and strengthening the Icelandic Low without weakening Aleutian.
  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v38 Southern Ocean wave checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Retuned the v30 shoulder-season warm-ocean source so its post-object
  amplitude depends more strongly on basin scale and open-ocean exposure.  This
  preserves a weak SON North Pacific source while reducing the v30 North
  Atlantic shoulder-season over-deepening.  MAM remains untriggered.
- Added a winter-only Northern Hemisphere subpolar source-expression boost.
  DJF Aleutian and Icelandic lows deepen modestly without changing shoulder
  seasons or Southern Hemisphere summer.
- Reworked the Southern Ocean gate: front, shelf/slope, and same-latitude SST
  anomaly now dominate source placement, while the longitude wave is only a
  weak perturbation.  Added a signed Southern Ocean wave-transfer anomaly so
  high-support sectors become relative lows and low-support sectors become
  relative highs.
- Rejected intermediate trials: v31 erased the SON North Pacific source; v32
  and v33 restored it while reducing North Atlantic excess; v35 and v36 were
  too conservative over the Southern Ocean; v37 introduced the correct signed
  wave-transfer direction; v38 makes that wave visible but bounded.
- This remains M2-only.  No R2b wind, R3 current, SST, precipitation, biome, or
  generated-world fitting was used.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v38_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- SON North Pacific keeps a weak local low source; SON North Atlantic is no
  longer over-deepened as strongly as v30.
- DJF Aleutian and Icelandic lows are modestly deeper than v30.
- JJA Southern Ocean now has a signed transfer wave.  The `60..120E` sector
  moves toward the real-Earth low, and several erroneous low sectors are
  lifted.  The wave is still not strong/organized enough for R2a promotion.
- Final pressure remains too smooth relative to Earth SLP, so v38 is a
  checkpoint, not an acceptance point.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.309/0.416/0.266`.
- Standardized-pressure correlation all/land/ocean:
  `0.564/0.674/0.272`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.588/0.600/0.576`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.11s`).
- Full test suite was started and interrupted for time after
  `103 passed in 444.62s`; no failure was observed before the interrupt.
Next:
- Stay in R2a/M2.  v38 is the current worktree checkpoint, but not an
  acceptance point.  The next owner is pressure-center compactness and
  remaining North Atlantic / Southern Ocean placement.  R2b/R3/R4/R5/R6/R7/R8
  and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v40 compact NH winter low checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Retained v38's shoulder-season source and Southern Ocean signed-wave logic.
- Strengthened only the Northern Hemisphere winter subpolar ocean-low
  compact-expression term, using existing low-support, SST-front, winter
  hemisphere, and subpolar-ocean gates.  This deepens existing
  Aleutian/Icelandic-type lows without activating MAM/SON/JJA ocean lows.
- v39 verified the direction at lower amplitude; offline replay showed that
  the stronger v40 coefficient improves North Pacific / North Atlantic winter
  core residuals and ocean pressure correlation before global MAE starts to
  flatten.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v40_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF Aleutian core replay deepens from `-0.575` in v38 to `-0.645` in v40
  against Earth `-0.840`.
- DJF Icelandic core deepens from `-0.435` in v38 to `-0.493` in v40 against
  Earth `-0.706`.
- MAM and SON North Pacific/North Atlantic source activation is unchanged by
  this winter-only term.  The v38 Southern Ocean signed transfer remains in
  place.
- R2a is still not promoted.  North Atlantic edge/Labrador placement and
  Southern Ocean wave organization remain incomplete, and final pressure is
  still smoother than Earth SLP.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.309/0.416/0.265`.
- Standardized-pressure correlation all/ocean:
  `0.565/0.277`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.589/0.577`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.26s`).
Next:
- Stay in R2a/M2.  v40 is now a historical checkpoint, not an
  acceptance point.  The next owner is North Atlantic edge/Labrador placement,
  remaining DJF pressure-center compactness, and Southern Ocean wave
  organization.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as
  fitting targets.

2026-07-07 - R2a M2 pressure-genesis v44 coastal inheritance and subpolar-front floor checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a narrow DJF coastal-land inheritance term for Northern Hemisphere
  subpolar ocean lows.  It lets Labrador/Iceland/coastal island land inherit
  adjacent ocean-low pressure without spreading into broad continental
  interiors.
- Added a bounded 50-67 N subpolar SST-front / storm-track low-support floor
  so Labrador/Icelandic winter lows are less suppressed by the earlier
  open-ocean exposure gate.
- Rejected v43: putting lee-low support into the pressure-object score
  overdeepened Labrador and NW Pacific while weakening the Icelandic low through
  object-selection competition.  v44 keeps the successful v42 coastal
  inheritance and uses a direct bounded floor instead.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v44_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF Labrador ocean residual improves from `+0.253` in v42 to `+0.179` in
  v44.
- DJF Icelandic ocean residual improves from `+0.353` in v42 to `+0.333` in
  v44.
- DJF Labrador land residual improves from `+0.242` in v42 to `+0.235` in
  v44; Icelandic land improves from `+0.512` to `+0.504`.
- NE Atlantic ocean becomes slightly overdeepened (`-0.014 -> -0.033`), so the
  next pass must avoid broad North Atlantic amplification.
- R2a is still not promoted.  The North Atlantic / Icelandic low remains too
  smooth and underdeep, Arctic/Nordic winter low structure is weak, and
  Southern Ocean sectors still need organization.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.308/0.415/0.265`.
- Standardized-pressure correlation all/ocean:
  `0.566/0.281`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.590/0.577`.
- DJF standardized-pressure MAE:
  `0.230`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.40s`).
Next:
- Stay in R2a/M2.  v44 is now a historical checkpoint, not an
  acceptance point.  The next owner is residual North Atlantic/Icelandic
  compactness and Southern Ocean wave-sector organization.  R2b/R3/R4/R5/R6/R7/R8
  and generated worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v45 Atlantic-Arctic gateway checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added an Atlantic-Arctic gateway low-support floor.  It computes graph
  distance from existing Atlantic subpolar low support through the
  Atlantic/Arctic ocean domain, then applies a bounded floor only to nearby
  62-80 N Arctic-basin marginal seas with shelf and SST-front support.
- This is a targeted gateway mechanism, not a blanket Arctic low: Pacific-side
  Arctic cells remain effectively unchanged.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v45_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF Nordic/Arctic gateway residual improves from `+0.642` in v44 to
  `+0.538` in v45.
- DJF Greenland Sea residual improves from `+0.645` to `+0.562`;
  Barents/Kara from `+0.479` to `+0.425`; Icelandic ocean from `+0.337` to
  `+0.314`.
- Beaufort/Arctic and Bering/Chukchi remain protected, so the fix does not
  create a spurious pan-Arctic low.
- R2a is still not promoted.  Nordic/Arctic and Icelandic lows remain too weak,
  and Southern Ocean shoulder-season / wave-sector organization is still
  visibly incomplete.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.308/0.415/0.265`.
- Standardized-pressure correlation all/ocean:
  `0.566/0.284`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.590/0.579`.
- DJF standardized-pressure MAE:
  `0.229`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.25s`).
Next:
- Stay in R2a/M2.  v45 is now a historical checkpoint, not an
  acceptance point.  The next owner is Southern Ocean shoulder-season /
  wave-sector source organization, with residual North Atlantic compactness as
  a secondary watch item.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain
  blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v46 Southern Ocean shoulder-source checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a Southern Ocean shoulder-season low-support floor.  It uses the
  existing Southern Ocean front/shelf/wave gate to create MAM/SON low-source
  support in the semantic Southern Ocean basin.
- This corrects the v45 failure mode where MAM/SON Southern Ocean had M1
  front/shelf support but zero M2 low-source expression.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v46_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- Semantic Southern Ocean MAM residuals improve in every 60-degree sector.
- Semantic Southern Ocean SON residuals improve in all target low-pressure
  sectors except that the 0-60E sector becomes more overdeepened.
- Wide 45-75S SON residuals improve at `-120..-60`, `60..120`, and
  `120..180`; `-60..60` remains too low.
- R2a is still not promoted.  v46 fixes the missing shoulder-source class but
  the Southern Ocean response is too band-like and needs wave-sector sharpening.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.306/0.416/0.262`.
- Standardized-pressure correlation all/ocean:
  `0.571/0.308`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.591/0.577`.
- MAM/SON standardized-pressure MAE:
  `0.355/0.430`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.58s`).
Next:
- Stay in R2a/M2.  v46 is the current worktree checkpoint, but not an
  acceptance point.  The next owner is Southern Ocean wave-sector sharpening,
  not broader Southern Ocean amplitude.  Residual North Atlantic compactness is
  a secondary watch item.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain
  blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v48 Southern Ocean shoulder-wave checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added bounded Southern Ocean shoulder-season transfer geometry on top of the
  v46 shoulder-source class.  The term uses semantic Southern Ocean support,
  latitude gates, SST-front support, shelf support, open-ocean exposure, and
  same-latitude SST anomaly.
- v47 added the first MAM/SON wave-sector terms.  It improved SON but pushed
  MAM too poleward.  v48 keeps the SON Pacific/Amundsen open-ocean low anchor
  and narrows the MAM transfer support toward the more equatorward storm-track
  belt.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v48_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- v48 makes the MAM/SON Southern Ocean transfer visibly more wave-sector based
  instead of a nearly uniform shoulder-season band.
- Wide 45-75S SON residuals improve from v46 in the Pacific/Amundsen,
  Atlantic-sector, Australian, and 120-180E sectors; the largest remaining
  positive residuals are still `-180..-120` and `120..180`.
- Semantic Southern Ocean SON residuals also improve in the target low-pressure
  sectors, but the -60..60 sector still has latitude-dependent sign conflict.
- MAM recovers from the v47 poleward overshoot and is slightly better than v46
  by MAM MAE and zonal-correlation.
- R2a is still not promoted.  The next owner remains Southern Ocean sector
  amplitude and latitude placement, with North Atlantic/Icelandic compactness
  as a secondary watch item.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.303/0.415/0.258`.
- Standardized-pressure correlation all/ocean:
  `0.580/0.355`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.599/0.597`.
- MAM/SON standardized-pressure MAE:
  `0.354/0.417`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.56s`).
Next:
- Stay in R2a/M2.  v48 is now the current worktree checkpoint, not an
  acceptance point.  Continue reading real-Earth SLP/replay/residual/source/
  transfer/M0/M1 maps before metrics.  R2b/R3/R4/R5/R6/R7/R8 and generated
  worlds remain blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v49 Southern Ocean latitude-split checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Split SON Southern Ocean wave-sector transfer into two latitude-dependent
  mechanisms: a subantarctic 45-62 S north-flank wave across southern
  Atlantic/Pacific/Indian margins, and a 58-78 S polar-side trough wave inside
  the semantic Southern Ocean.
- The split fixes the v48 problem where one phase had to represent both the
  north-flank ridge and the polar-side trough.  Both terms remain bounded and
  tied to SST-front, shelf/open-ocean, same-latitude SST, and semantic basin
  support.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v49_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- SON Southern Ocean transfer now has readable latitude-dependent wave
  geometry rather than a single annular shoulder response.
- Wide 45-75S SON residuals improve in nearly all sectors, especially
  `-180..-120`, `-60..60`, and `120..180`.
- Semantic Southern Ocean `lat < -55` residuals improve strongly at
  `120..180` and moderately at `-180..-120`.
- R2a is still not promoted.  Remaining visible owners are high-latitude
  Southern Ocean `-180..0` trough underdepth, DJF North Atlantic / Nordic /
  Barents gateway compactness, and muted JJA subpolar-ocean pressure phase.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.298/0.412/0.252`.
- Standardized-pressure correlation all/ocean:
  `0.593/0.402`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.612/0.625`.
- MAM/SON standardized-pressure MAE:
  `0.354/0.395`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.13s`).
Next:
- Stay in R2a/M2.  v49 is now the current worktree checkpoint, not an
  acceptance point.  Select the next owner by real-Earth map read; do not
  resume R2b/R3/R4/R5/R6/R7/R8 or generated-world fitting yet.

2026-07-07 - R2a M2 pressure-genesis v50 Atlantic-Arctic gateway transfer checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a bounded DJF Atlantic-Arctic gateway transfer term for the Iceland /
  Greenland / Norwegian / Barents gateway.  The term uses latitude, longitude,
  Atlantic/Arctic basin, shelf, and SST-front support.
- Added a small coastal-land inheritance term from the same gateway seed so
  Iceland/Greenland/Nordic coastal cells can inherit nearby low pressure.
- The term avoids Beaufort and Bering/Chukchi longitudes, preserving protection
  against a spurious pan-Arctic low.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v50_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- DJF Greenland Sea, Norwegian Sea, Barents/Kara, Icelandic core, and Iceland
  coastal residuals all move toward Earth.
- Beaufort and Bering/Chukchi remain nearly unchanged, so v50 does not create a
  blanket Arctic low.
- R2a is still not promoted.  The next visible owner is warm-season /
  shoulder-season subpolar ocean positive-pressure support: current M2
  high-support objects are too confined to the subtropics.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.297/0.412/0.251`.
- Standardized-pressure correlation all/ocean:
  `0.594/0.406`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.613/0.627`.
- DJF standardized-pressure MAE:
  `0.227`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.84s`).
Next:
- Stay in R2a/M2.  v50 is now the current worktree checkpoint, not an
  acceptance point.  The next map-read owner should be MAM/JJA subpolar ocean
  high-pressure support, while preserving Southern Ocean and gateway fixes.

2026-07-07 - R2a M2 pressure-genesis v51 subpolar ocean high-support checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added MAM/JJA Northern Hemisphere subpolar-ocean high-pressure support.
  MAM support is broad over high-latitude ocean/shelf/frontal seas; JJA support
  is narrower and longitude-gated to North Pacific / Gulf of Alaska / North
  Atlantic margins while protecting Beaufort.
- The support is exported through `ocean_pressure_high_source_support` and
  receives bounded positive M2 source-to-pressure transfer.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v51_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/real_earth_pressure_replay_contact_sheet.png`.
- High-support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/replay_ocean_pressure_high_source_support_seasons.png`.
Map read:
- The missing MAM/JJA subpolar ocean high object class is now present.
- MAM North Pacific / North Atlantic / Barents / Bering residuals move strongly
  toward Earth.  Greenland Sea and Beaufort remain under-high but no longer
  lack an object class.
- JJA North Pacific, Gulf of Alaska, North Atlantic, Labrador, and Bering
  improve, while Beaufort remains protected.
- R2a is still not promoted.  Remaining owners include MAM polar-cap
  regularity and residual Greenland/Beaufort under-high, JJA North Pacific /
  Gulf of Alaska under-high, and land shoulder-season pressure errors.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.290/0.411/0.240`.
- Standardized-pressure correlation all/ocean:
  `0.614/0.499`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.611/0.621`.
- MAM/JJA standardized-pressure MAE:
  `0.329/0.207`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.71s`).
Next:
- Stay in R2a/M2.  v51 is now the current worktree checkpoint, not an
  acceptance point.  Select the next owner by real-Earth map read; do not
  resume downstream wind/current/SST/precipitation/biome or generated-world
  fitting.

2026-07-07 - R2a M2 pressure-genesis v52 land shoulder-phase checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a Northern Hemisphere high-latitude land shoulder-season phase
  correction inside M2 pressure genesis.  MAM receives a residual cold/snow-
  shield high-pressure support from spring cold anomaly, near-freezing
  temperature, low interiority, coast strength, and latitude gates.  SON
  receives a warm-ground / summer-heat-memory high-pressure decay support that
  suppresses the previous premature autumn continental high.
- The correction is expressed through both land pressure source support and
  bounded M2 source-to-pressure transfer.  It does not change the upstream
  geographic circulation / wind / ocean-current generator.
- R2b wind, R3 currents, SST, precipitation, biomes, generated worlds, and
  downstream acceptance metrics remain blocked as fitting targets.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v52_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Land support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_land_pressure_source_support_seasons.png`.
- Transfer map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_pressure_genesis_wave_transfer_seasons.png`.
Map read:
- MAM land shoulder residuals move strongly toward Earth: Siberia
  `-0.487 -> -0.114`, Eurasia high-latitude land `-0.442 -> +0.025`,
  Canada `-1.145 -> -0.642`, North America high-latitude land
  `-1.013 -> -0.469`, and Greenland `-0.715 -> +0.047`.
- SON premature land high residuals also move strongly toward Earth:
  Siberia `+1.250 -> +0.328`, Eurasia high-latitude land `+0.909 -> -0.046`,
  Canada `+1.147 -> +0.109`, North America high-latitude land
  `+1.072 -> +0.177`, and Greenland `+0.057 -> -0.033`.
- The MAM polar-cap / Arctic-edge red band remains visible and was already
  present in v51; v52 does not solve that old blocker.
- R2a is still not promoted.  Remaining owners are MAM polar-cap
  latitude/texture regularity and Greenland/Beaufort under-high, JJA North
  Pacific / Gulf of Alaska under-high, and residual Southern Ocean high-lat
  trough underdepth.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.284/0.394/0.239`.
- Standardized-pressure correlation all/land/ocean:
  `0.626/0.689/0.495`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.597/0.590/0.612`.
- MAM/SON standardized-pressure MAE:
  `0.313/0.388`.
- MAM/SON zonal-anomaly correlation:
  `0.566/0.436`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.73s`).
Next:
- Stay in R2a/M2.  v52 is the current worktree checkpoint, not an acceptance
  point.  The next map-read owner should be MAM polar-cap regularity /
  Greenland-Beaufort under-high or JJA North Pacific / Gulf of Alaska under-
  high.  Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic
  gateway, v51 Beaufort-protected subpolar-ocean high support, and v52 land
  shoulder-phase correction.

2026-07-07 - R2a M2 pressure-genesis v53 MAM Arctic freeze-ocean high checkpoint, not promoted
Changed:
- Kept the work inside the real-Earth R2a seasonal SLP / pressure-source
  packet.
- Added a MAM Arctic / Greenland / Beaufort freeze-ocean high-pressure object
  in M2 pressure genesis.  The support is derived from near-freezing SST,
  Arctic/Baffin/Greenland/Beaufort longitude gates, shelf support, SST-front
  support, and high-latitude ocean masks.
- Added a Baffin/Labrador gateway sub-support because that region is lower
  latitude and was under-triggered by the polar cap gate alone.
- The mechanism contributes to both ocean high support/source and bounded
  source-to-pressure transfer.  It does not alter wind/current/SST generation.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v53_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/real_earth_pressure_replay_contact_sheet.png`.
- High-support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/replay_ocean_pressure_high_source_support_seasons.png`.
Map read:
- MAM Beaufort residual improves `-0.893 -> -0.242`.
- MAM Greenland Sea residual improves `-0.661 -> -0.035`.
- MAM Barents-Kara residual changes `-0.239 -> +0.105`, acceptable for this
  checkpoint.
- MAM Baffin/Labrador residual improves `-0.532 -> -0.167`.
- MAM Arctic cap residual improves `-0.662 -> -0.262`.
- R2a remains unpromoted.  Canada MAM land remains too low, JJA North Pacific
  / Gulf Alaska remains under-high, and SON remains the weakest season.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.281/0.391/0.236`.
- Standardized-pressure correlation all/land/ocean:
  `0.631/0.684/0.522`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.601/0.591/0.620`.
- MAM standardized-pressure MAE / zonal-anomaly correlation:
  `0.300/0.587`.

2026-07-07 - R2a M2 pressure-genesis v54 JJA North Pacific high checkpoint, not promoted
Changed:
- Added a JJA North Pacific / Gulf of Alaska high-pressure object in M2
  pressure genesis.  The support is restricted to semantic Pacific basin
  ocean, 41-67 N, and Gulf Alaska / Aleutian longitude gates, with shelf,
  SST-front, and cool same-latitude SST support.
- The support is intentionally zero over Beaufort/Arctic so v51's
  Beaufort-protection remains intact.
- This remains R2a/M2-only; R2b wind, R3 currents, SST, precipitation, biomes,
  and generated worlds are still not fitting targets.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v54_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/real_earth_pressure_replay_contact_sheet.png`.
- High-support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/replay_ocean_pressure_high_source_support_seasons.png`.
Map read:
- JJA Gulf Alaska residual improves `-0.585 -> -0.268`.
- JJA Aleutian residual improves `-0.446 -> -0.156`.
- JJA North Pacific residual improves `-0.376 -> -0.121`.
- JJA NW Pacific changes `-0.117 -> +0.139`; this is a mild over-high tradeoff
  but still within the intended North Pacific high-pressure object.
- JJA Beaufort remains protected (`+0.308 -> +0.303`).
- R2a remains unpromoted.  Remaining visible owners are MAM Canada land
  under-high, residual Arctic/Beaufort MAM under-high, SON pressure wave /
  Southern Ocean residuals, and overall pressure-center texture.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.279/0.391/0.234`.
- Standardized-pressure correlation all/land/ocean:
  `0.633/0.684/0.530`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.604/0.593/0.624`.
- MAM standardized-pressure MAE / zonal-anomaly correlation:
  `0.300/0.587`.
- JJA standardized-pressure MAE / zonal-anomaly correlation:
  `0.203/0.741`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.66s`).
Next:
- Stay in R2a/M2.  v54 is the current worktree checkpoint, not an acceptance
  point.  Select the next owner by real-Earth map read, likely MAM Canada land
  under-high or SON / Southern Ocean pressure-wave residuals.  Continue to
  preserve v49-v54 pressure-source fixes and keep downstream systems blocked.

2026-07-07 - R2a M2 pressure-genesis v55 MAM North America land-high checkpoint, not promoted
Changed:
- Added a MAM North America spring land-high object in M2 pressure genesis.
  The support is restricted to 43-73 N and -165..-45 longitude, with a compact
  Canada-centered longitude gate, low-elevation support, cold / near-freezing
  memory, low-interiority weighting, and coast-strength modulation.
- This is intentionally narrower than the generic land shoulder-phase term:
  the generic candidate also raised Eurasia/Siberia, while this object targets
  the remaining real-Earth Canada / North America under-high.
- Fixed the land-support diagnostic merge to use `np.maximum.reduce`, removing
  the previous NumPy warning from a three-argument `np.maximum` call.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v55_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Land-support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/replay_land_pressure_source_support_seasons.png`.
Map read:
- MAM Canada residual improves `-0.665 -> -0.321`.
- MAM North America high-latitude land residual improves `-0.502 -> -0.243`.
- MAM Greenland land remains near target (`-0.054 -> -0.009`).
- MAM Siberia / Eurasia remain near target (`-0.107 -> -0.123`,
  `+0.008 -> -0.009`), so the fix does not reintroduce broad Eurasian
  over-high.
- MAM Alaska becomes mildly over-high (`+0.066 -> +0.244`), retained as an
  acceptable tradeoff because the broader Canada/North America land object and
  MAM metrics improve.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.278/0.387/0.234`.
- Standardized-pressure correlation all/land/ocean:
  `0.635/0.688/0.529`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.607/0.600/0.623`.
- MAM standardized-pressure MAE / zonal-anomaly correlation:
  `0.296/0.604`.
Tests:
- Targeted regression tests pass without the previous NumPy warning:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.57s`).
Next:
- Stay in R2a/M2.  v55 is the current worktree checkpoint, not an acceptance
  point.  The next visible owner is likely SON / Southern Ocean pressure-wave
  residuals or residual MAM Arctic cap under-high.  Preserve v49-v55
  pressure-source fixes and keep downstream systems blocked.

2026-07-07 - R2a M2 pressure-genesis v58 DJF ocean lows and North America winter-high relief checkpoint, not promoted
Changed:
- Added DJF Atlantic-Arctic gateway low support and DJF North Pacific /
  Aleutian low support in M2 pressure genesis.  These terms are gated by
  semantic basin IDs, latitude/longitude, shelf support, and SST-front support,
  and are expressed through bounded negative source plus source-to-pressure
  transfer.  They do not alter R2b wind, R3 currents, SST, precipitation, or
  generated-world fitting targets.
- Strengthened the Atlantic-Arctic branch after v56 showed Icelandic/Nordic
  improvement but residual Greenland Sea / Barents under-depth.
- Added DJF North America winter-high relief support: a bounded negative
  pressure-source object over North America, representing weaker continental
  winter high pressure than Siberia because the smaller continent is more
  exposed to Pacific/Atlantic maritime erosion and adjacent storm-track lows.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v58_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Residual map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/pressure_standardized_delta_seasons.png`.
- Land-support map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/replay_land_pressure_source_support_seasons.png`.
Map read:
- DJF Icelandic/Nordic residual improves `+0.188 -> +0.043`.
- DJF Greenland Sea residual improves `+0.400 -> +0.145`; still not perfect,
  but no longer the dominant winter low failure.
- DJF Barents/Kara improves `+0.216 -> +0.051`; Labrador/Baffin improves
  `+0.158 -> +0.038`.
- DJF Aleutian improves `+0.112 -> +0.050`; Bering/Chukchi remains bounded
  (`-0.048 -> -0.092`) rather than over-deepened.
- DJF Canada land over-high improves `+0.425 -> +0.244`, while Siberia remains
  near target (`-0.001 -> +0.023`).
- R2a remains unpromoted.  Visible remaining owners are MAM Canada land
  under-high, MAM Arctic cap under-high, JJA Gulf Alaska under-high, and the
  broader SON residual pattern.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.277/0.385/0.234`.
- Standardized-pressure correlation all/land/ocean:
  `0.638/0.690/0.534`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.610/0.603/0.626`.
- DJF standardized-pressure MAE / zonal-anomaly correlation:
  `0.222/0.720`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.94s`).
Next:
- Stay in R2a/M2.  v58 is the current worktree checkpoint, not an acceptance
  point.  Select the next owner by real-Earth map read, likely MAM Arctic /
  Canada seasonal high-pressure support or JJA Gulf Alaska high support.  Keep
  R2b wind, R3 currents, SST, precipitation, biomes, and generated worlds
  blocked as fitting targets.

2026-07-07 - R2a M2 pressure-genesis v59 MAM Canada / Canadian Arctic high checkpoint, not promoted
Changed:
- Added a MAM central/eastern North America plains high object.  The support
  is restricted away from Alaska/Yukon and gated by latitude/longitude,
  low-elevation plains/shield terrain, low interiority, coast-strength memory,
  and terrain shelter.
- Added MAM Canadian Arctic and central-west Arctic freeze-high objects.  These
  use Arctic semantic basin support, near-freezing SST, latitude/longitude
  gates, shelf / deep-basin weighting, and explicit Barents/Kara protection.
- The changes are M2 pressure-source geometry only.  R2b wind, R3 currents,
  SST, precipitation, biomes, and generated worlds remain blocked.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v59_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Residual map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/pressure_standardized_delta_seasons.png`.
Map read:
- MAM Canada land residual improves `-0.309 -> -0.015`.
- MAM west / central / east Canada improve `-0.223 -> -0.012`,
  `-0.570 -> -0.261`, and `-0.275 -> +0.009`.
- MAM lower North America improves `-0.298 -> +0.030`.
- MAM Alaska/Yukon is not worsened (`+0.258 -> +0.226`), so the new plains
  object avoids the prior over-high tradeoff.
- MAM Arctic cap improves `-0.289 -> -0.231`; Beaufort improves
  `-0.266 -> -0.068`; Canadian Archipelago improves `-0.475 -> -0.255`.
- MAM Barents/Kara improves `+0.097 -> +0.058`; the new Arctic support did
  not create a Barents over-high wall.  Baffin/Labrador and central Arctic
  remain under-high and are possible later owners.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.276/0.381/0.233`.
- Standardized-pressure correlation all/land/ocean:
  `0.640/0.691/0.539`.
- Pressure zonal-anomaly correlation all/land/ocean:
  `0.613/0.607/0.626`.
- MAM standardized-pressure MAE / zonal-anomaly correlation:
  `0.290/0.615`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.31s`).
Next:
- Stay in R2a/M2.  v59 is the current worktree checkpoint, not an acceptance
  point.  Next owner should be chosen by map read between residual
  Baffin/Labrador + central Arctic MAM under-high, JJA Gulf Alaska / Aleutian
  under-high, and broader SON pressure-wave residuals.

2026-07-07 - R2a M2 pressure-genesis v60 JJA eastern / central North Pacific high checkpoint, not promoted
Changed:
- Added a JJA eastern / central North Pacific high-pressure object in M2
  pressure genesis.  The support is restricted to semantic Pacific basin cells
  in the Gulf Alaska / eastern Aleutian / central North Pacific sector and is
  gated by latitude, longitude, shelf/front support, cool same-latitude SST,
  and an Arctic taper.
- The object explicitly suppresses the already over-high NW Pacific
  130E-160E sector and does not touch Beaufort/Arctic or North Atlantic
  pressure support.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v60_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Residual map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/pressure_standardized_delta_seasons.png`.
Map read:
- JJA Gulf Alaska residual improves `-0.265 -> -0.118`.
- JJA Gulf Alaska east / west improve `-0.235 -> -0.089` and
  `-0.309 -> -0.170`.
- JJA eastern Aleutian improves `-0.343 -> -0.181`; western Aleutian changes
  `-0.101 -> +0.014`, which is an acceptable local tradeoff.
- JJA North Pacific 160W-130W and 180-160W improve `-0.247 -> -0.118` and
  `-0.251 -> -0.136`.
- NW Pacific 130E-160E remains over-high but does not worsen
  (`+0.369 -> +0.360`); Beaufort remains protected (`+0.303 -> +0.299`).
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.275/0.381/0.232`.
- Standardized-pressure correlation all/ocean:
  `0.642/0.543`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.615/0.631`.
- JJA standardized-pressure MAE / zonal-anomaly correlation:
  `0.199/0.749`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.41s`).
Next:
- Stay in R2a/M2.  v60 is the current worktree checkpoint, not an acceptance
  point.  Next visible owners are residual MAM Baffin-Labrador / central
  Arctic under-high, JJA North Atlantic subpolar under-high, and the broader
  SON pressure-wave residual pattern.

2026-07-07 - R2a M2 pressure-genesis v62 SON shoulder-season pressure checkpoint, not promoted
Changed:
- Added a SON boreal autumn land-relief support in M2 pressure genesis.  It
  damps early / over-broad low-elevation continental cold highs using
  summer-heat memory, autumn cooling, unfrozen-ground memory, low-elevation
  support, coast / low-interiority leakage, and terrain-barrier escape.  This
  targets the v60 SON Siberia / Canada / Alaska over-high without using wind,
  ocean-current, SST, precipitation, biome, or generated-world fitting.
- Added a SON North Atlantic autumn low support.  It is restricted to semantic
  Atlantic ocean cells around the Icelandic-low sector and gated by latitude,
  longitude, shelf, SST-front, and same-latitude SST support; a west taper
  avoids over-deepening Labrador.
- Added a small SON Southern Ocean sector low wave adjustment after map-read
  showed the v61 land / North Atlantic fix raised Southern Ocean residuals.
  The adjustment is strongest in Pac/Amundsen and Indian sectors and uses an
  Atlantic cut so the Atlantic sector stays near target.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v62_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Residual map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/pressure_standardized_delta_seasons.png`.
Map read:
- SON Siberia residual improves `+0.338 -> +0.137`; Canada improves
  `+0.204 -> -0.057`; Alaska improves `+0.395 -> +0.078`.
- SON North Atlantic subpolar improves `+0.247 -> +0.074`; Icelandic sector
  improves `+0.455 -> +0.098`; Labrador is preserved (`+0.076 -> +0.078`).
- SON Southern Ocean all returns to the v60 level (`+0.066 -> +0.066`) after
  the v61 side effect.  Pac/Amundsen is near v60 (`+0.106 -> +0.109`);
  Indian improves relative to v61 (`+0.211 -> +0.133`) and improves relative
  to v60 (`+0.183 -> +0.133`).
- MAM central Arctic / Baffin-Labrador and JJA North Atlantic / NW Pacific are
  intentionally unchanged; they remain later R2a owners.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.272/0.377/0.230`.
- Standardized-pressure correlation all/land/ocean:
  `0.644/0.693/0.547`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.611/0.624`.
- SON standardized-pressure MAE / zonal-anomaly correlation:
  `0.378/0.416`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.30s`).
Next:
- Stay in R2a/M2.  v62 is the current worktree checkpoint, not an acceptance
  point.  Remaining visible owners are MAM Baffin-Labrador / central Arctic
  under-high, JJA North Atlantic subpolar under-high and NW Pacific over-high,
  plus residual SON high-latitude texture / Antarctica edge artifacts.

2026-07-07 - R2a M2 pressure-genesis v63 JJA Eurasian thermal-low / western Pacific trough checkpoint, not promoted
Changed:
- Added a JJA Eurasian summer thermal-low object in M2 pressure genesis.  It
  is restricted to the heated low-elevation Eurasian land belt from Arabia /
  Iran through India and East Asia and is gated by JJA thermal anomaly,
  low-elevation support, continent / coast expression, and terrain-barrier
  escape.  This addresses an underexpressed continental summer low without
  touching North America, which was already too low.
- Added a JJA western Pacific marginal-sea low / trough object.  It is
  restricted to semantic Pacific / Indian marginal-sea cells in the 118E-178E,
  24-64N domain and is gated by shelf, SST front, and same-latitude SST
  support.  This targets the high-pressure sign error over the Kuroshio /
  Oyashio and Japan / East China Sea sector while preserving the Gulf Alaska
  and central North Pacific high.
Evidence:
- Real-Earth replay output:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v63_20260707/`.
- Pressure replay output:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/`.
- Main contact sheet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/real_earth_pressure_replay_contact_sheet.png`.
- Residual map:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/pressure_standardized_delta_seasons.png`.
Map read:
- JJA NW Pacific residual improves `+0.360 -> +0.070`.
- JJA Kuroshio / Oyashio improves `+0.312 -> +0.065`.
- JJA Japan / East China Sea improves `+0.747 -> +0.606`, but remains a
  visible residual owner.
- JJA East Asia land improves `+0.386 -> +0.225`; NE Asia improves
  `+0.223 -> +0.072`; China lowland improves `+0.539 -> +0.352`.
- JJA India improves `+0.431 -> +0.312`; Arabia / Iran improves
  `+0.636 -> +0.499` but remains under-low.
- JJA Gulf Alaska is preserved (`-0.132 -> -0.128`), central North Pacific is
  preserved (`-0.172 -> -0.153`), and MAM / SON target regions are unchanged.
Metrics:
- Standardized-pressure MAE all/land/ocean:
  `0.271/0.373/0.229`.
- Standardized-pressure correlation all/land/ocean:
  `0.649/0.699/0.551`.
- Pressure zonal-anomaly correlation all/ocean:
  `0.618/0.630`.
- JJA standardized-pressure MAE / zonal-anomaly correlation:
  `0.192/0.778`.
Tests:
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 105.55s`).
Next:
- Stay in R2a/M2.  v63 is the current worktree checkpoint, not an acceptance
  point.  Remaining visible owners are MAM Baffin-Labrador / central Arctic
  under-high, JJA North Atlantic / Icelandic under-high, residual JJA Japan /
  East China Sea and Arabia / Iran summer low underexpression, and residual
  SON high-latitude texture / Antarctica edge artifacts.
