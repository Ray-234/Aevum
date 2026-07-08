# R2a M0/M1 Map-Read Attribution

Status: active R2a attribution note
Date: 2026-07-07

This note records the first real-Earth visual comparison after adding the M0/M1
diagnostic contract.  It follows the Earth-only replay protocol: compare the
real-Earth subgraph with the Aevum same-grid replay first, then assign visible
residuals to an owner before changing code.

## Evidence Packet

Earth reference:

- `out_earth_climate_reference_r6_landcover_20260706/earth_reference_8000cells.npz`

Aevum replay:

- `out_real_earth_climate_replay_m0_m1_contract_20260707/`
- `out_real_earth_climate_replay_m0_m1_contract_20260707/terminal_climate_arrays.npz`

Visual comparison:

- `out_real_earth_pressure_replay_m0_m1_contract_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_m0_m1_contract_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_m0_m1_contract_20260707/replay_m1_energy_support_contact_sheet.png`

Regression metrics:

- pressure standardized MAE all/land/ocean:
  `0.317 / 0.417 / 0.277`
- pressure standardized corr all/land/ocean:
  `0.539 / 0.668 / 0.135`
- pressure zonal-anomaly corr all/land/ocean:
  `0.567 / 0.578 / 0.560`
- wind speed MAE:
  `1.865 m/s`
- wind direction cosine p50:
  `0.858`

The pressure metrics match the previous R2a diagnostic baseline.  This is
expected: the M0/M1 implementation is diagnostic-only and does not change the
accepted climate physics.

## Visual Read

Real Earth SLP:

- Earth has basin-scale winter lows and highs, especially North Pacific /
  Aleutian, North Atlantic / Icelandic, and seasonally shifted continental
  highs/lows.
- The Southern Ocean pressure field has a strong annular component, but it is
  not only a single featureless stripe.  It carries wave-like and basin-adjacent
  structure.
- Coast, mountain, and land-sea stationary-wave structures are visible in the
  seasonal residual patterns.

Aevum pressure replay:

- The replay captures broad continental winter highs and summer lows.
- It still reads as large smooth continental blobs plus a strong Southern Ocean
  latitude band.
- Ocean pressure-center structure is weak.  North Pacific / North Atlantic
  ocean pressure centers do not organize like Earth.
- Stationary-wave support is visually present, but pressure centers do not yet
  convert it into Earthlike pressure geometry.

M0 boundary support:

- Elevation, barrier belts, shelf/coast, and lapse-cooling support maps are
  readable and largely useful for attribution.
- The major upstream defect is `ocean.basin_id`: in this Earth replay it is
  dominated by one connected global ocean basin.  Basin 0 contains about
  `99.2%` of ocean area, and the remaining ids are tiny enclosed/separated
  fragments.
- Because Atlantic, Pacific, Indian, Arctic, and Southern Ocean basin ids are
  not represented as major semantic basins, any R2a pressure mechanism that
  depends on `ocean.basin_id` cannot create realistic basin-specific pressure
  centers.

M1 energy support:

- Seasonal insolation and land thermal anomalies are present and physically
  interpretable.
- Land thermal anomalies and land-sea thermal contrast support continental
  winter highs and summer lows.  Therefore the broad continental pressure
  signal is not purely missing upstream energy.
- Ocean mixed-layer seasonal anomaly remains strongly latitude-controlled.
- Same-latitude SST anomaly and SST-gradient support do contain useful
  non-zonal ocean/coastal structure, but current M2 pressure genesis does not
  make enough basin-scale pressure-center use of them.
- Snow/ice albedo support is concentrated at high latitudes and can explain
  some polar residuals, but R7 feedback remains blocked until R2a-R6 pass.

## Residual Ownership

| Residual | Owner | Reason |
| --- | --- | --- |
| North Pacific / North Atlantic winter lows missing or weak | M0 then M2 | M0 lacks major ocean-basin semantic ids; M2 also does not yet consume SST-front support as pressure-center forcing. |
| Southern Ocean over-zonal pressure residual | M2, with M0 basin semantics as prerequisite | Current pressure source is too latitude-band based; basin/front perturbations are too weak. |
| Continental pressure blobs too smooth/blocky | M2 | M1 land thermal contrast exists, but pressure genesis smooths it into broad continent-scale patches. |
| Mountain/coast stationary-wave support present but under-expressed | M2 | M0 barrier/coast support exists; pressure-center geometry does not yet anchor enough to those supports. |
| R2b wind residuals | blocked by R2a | Wind must not be tuned while pressure map remains visibly wrong. |

## Acceptance

R2a is not accepted.

The next code change should not tune R2b wind, ocean currents, SST closure,
moisture, precipitation, ice, Koppen, or biomes.  It should repair the active
R2a owner path:

1. Add major-ocean semantic basin support for Earth replay without depending on
   accidental connected components.
2. Re-render M0 boundary support and verify Atlantic/Pacific/Indian/Arctic /
   Southern basin readability.
3. Add or adjust M2 pressure genesis so ocean basin/front support can create
   basin-scale Aleutian/Icelandic/Southern Ocean pressure centers.
4. Re-run the same pressure contact sheet and compare visually before reading
   metrics.

## Major-Ocean Basin Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_major_ocean_basins_20260707/`
- `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/replay_m1_energy_support_contact_sheet.png`

M0 basin-id result:

| Basin id | Semantic basin | Ocean area share | Cell count |
| --- | --- | ---: | ---: |
| 0 | Atlantic Ocean | `0.2477` | `1412` |
| 1 | Pacific Ocean | `0.3887` | `2216` |
| 2 | Indian Ocean | `0.2354` | `1342` |
| 3 | Arctic Ocean | `0.0375` | `214` |
| 4 | Southern Ocean | `0.0907` | `517` |

Map read:

- The M0 support sheet now has readable Atlantic, Pacific, Indian, Arctic, and
  Southern Ocean semantic ids.
- This resolves the specific upstream M0 defect where the Earth replay ocean
  collapsed into one connected component.
- The R2a pressure contact sheet is still not acceptable.  Replay pressure is
  still too smooth over continents, and the Southern Ocean still reads too much
  like a broad zonal band.  North Pacific and North Atlantic pressure-center
  geometry remains under-expressed.

Updated residual ownership:

| Residual | Previous owner | Current owner | Reason |
| --- | --- | --- | --- |
| North Pacific / North Atlantic winter lows missing or weak | M0 then M2 | M2 | Major ocean ids are now readable; pressure genesis still does not construct basin-scale low objects from basin/SST-front support. |
| Southern Ocean over-zonal pressure residual | M2, with M0 prerequisite | M2 | Southern Ocean semantic support exists, but pressure genesis still overuses latitude-band forcing. |
| Continental pressure blobs too smooth/blocky | M2 | M2 | M1 thermal support exists; pressure genesis smooths it into broad blobs. |
| Mountain/coast stationary-wave support under-expressed | M2 | M2 | M0 barrier/coast support exists; pressure objects do not yet anchor enough to it. |
| R2b wind residuals | blocked by R2a | blocked by R2a | Wind must remain frozen until R2a map-read acceptance. |

Next code step:

1. Keep the real-Earth major-ocean M0 support.
2. Build M2 pressure-center candidates from `ocean.basin_id`,
   `climate.same_latitude_sst_anomaly`, `climate.sst_gradient_support`,
   `atmosphere.stationary_wave_pressure_support`, `terrain.barrier_index`, and
   seasonal hemisphere forcing.
3. Render Earth SLP, Aevum pressure replay, residual, pressure-center support,
   and M0/M1 support again.
4. Judge by map read first: Aleutian Low, Icelandic Low, Southern Ocean wave
   structure, and continental winter-high/summer-low geometry must become more
   Earthlike before any R2b/R3/R4 work resumes.

## M2 Pressure-Genesis v1 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v1_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v1_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v1_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v1_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation:

- Added an object-based M2 pressure-source pass in `ClimateModule`.
- The pass extracts winter subpolar ocean-low candidates from major ocean basin
  ids, open-ocean support, SST-front support, and same-latitude SST anomaly.
- It adds weak basin subtropical highs and small land/terrain stationary-wave
  refinements.
- It updates `atmosphere.land_sea_pressure_proxy`; it does not retune R2b wind
  translation.

Map read:

- Improvement: DJF North Pacific and North Atlantic low-pressure objects are
  more visible than in the M0-only checkpoint.
- Improvement: JJA Southern Ocean pressure support is less like one completely
  uniform latitude band and starts to segment into basin/front-supported
  patches.
- Remaining defect: continents still read as overly smooth broad pressure
  blocks, especially in shoulder seasons.
- Remaining defect: Southern Ocean remains too zonal, and pressure-center
  amplitudes/positions are still not close enough to Earth SLP for promotion.
- Remaining defect: pressure-center support is useful but still partly reflects
  the broad base pressure field rather than clean Aleutian/Icelandic/Southern
  Ocean objects.

Metrics versus the major-ocean M0 checkpoint:

| Metric | M0 basin checkpoint | M2 v1 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.320` |
| standardized-pressure corr all | `0.539` | `0.545` |
| standardized-pressure corr ocean | `0.137` | `0.213` |
| pressure zonal-anomaly corr all | `0.568` | `0.577` |
| pressure zonal-anomaly corr ocean | `0.563` | `0.569` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |

Acceptance:

R2a is still not accepted.

Next M2 step:

1. Keep object extraction, but refine source placement and amplitude so ocean
   lows improve map geometry without increasing MAE as much.
2. Separate pressure-center support into causal source support versus residual
   broad pressure support, so map reading can distinguish generated Aleutian /
   Icelandic / Southern Ocean objects from inherited base blobs.
3. Improve continental stationary waves by using terrain orientation and
   continent-edge geometry, not just scalar land thermal anomaly.
4. Re-render the same Earth-only pressure packet before touching R2b wind.

## M2 Pressure-Genesis v2 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v2_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v1:

- The signed M2 pressure-source increment is archived as
  `atmosphere.pressure_genesis_source`.
- Ocean low, ocean high, land thermal, and terrain-wave causal supports are
  archived separately.
- Ocean-low source extraction is narrower and lower-amplitude than v1, reducing
  the MAE penalty while preserving basin/front source readability.

Map read:

- DJF North Pacific and North Atlantic ocean-low causal sources are directly
  visible in the new source panels.
- JJA Southern Ocean source support is visible and tied to the Southern Ocean
  basin/front band, but it is still too continuous and not wave-like enough.
- Final replay pressure is still smoother than Earth SLP, especially over large
  continents and shoulder seasons.
- Therefore R2a remains unaccepted.  The pressure source is now auditable, but
  source placement and continental/terrain stationary waves still need M2 work.

Metrics versus the M0 basin checkpoint:

| Metric | M0 basin checkpoint | M2 v2 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE ocean | `0.276` | `0.277` |
| standardized-pressure corr all | `0.539` | `0.546` |
| standardized-pressure corr ocean | `0.137` | `0.205` |
| pressure zonal-anomaly corr all | `0.568` | `0.576` |
| pressure zonal-anomaly corr ocean | `0.563` | `0.570` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |

Next M2 step:

1. Keep causal source diagnostics.
2. Break the Southern Ocean source into wave-like sectors using SST-front
   phase, terrain/coast downstream anchors, or basin-sector phase rather than
   accepting a continuous annular band.
3. Add continent-edge and terrain-orientation pressure-wave placement so
   continental highs/lows become less blob-like.
4. Continue to render the same Earth-only pressure packet; R2b/R3/R4 remain
   blocked.

## M2 Pressure-Genesis v4 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v4_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v2:

- Southern Ocean ocean-low extraction uses 60-degree sector labels plus a
  wavenumber-3/front gate.  The gate is deterministic and derived from SST
  front support plus reduced stationary-wave phase, not random blobs.
- Continental pressure refinement now uses continent-level thermal-center
  objects.  Large continents can keep multiple centers; terrain/land-thermal
  gradient support replaces the previous broad all-land scalar increment.
- R2b wind, R3 currents, SST, precipitation, and biomes remain unchanged.

Map read:

- DJF North Pacific and North Atlantic causal low sources remain readable.
- JJA Southern Ocean source is now visibly sectorized in the source panel,
  rather than a single uniform annular band.
- Continental source maps show object-like interior/terrain refinements, but
  the final pressure replay remains too smooth compared with real Earth SLP.
- R2a remains unaccepted.  The source layer is now much more auditable, but the
  final pressure field still needs stronger stationary-wave expression before
  R2b wind can be unfrozen.

Metrics versus the M0 basin checkpoint:

| Metric | M0 basin checkpoint | M2 v4 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE ocean | `0.276` | `0.276` |
| standardized-pressure corr all | `0.539` | `0.547` |
| standardized-pressure corr land | `0.668` | `0.668` |
| standardized-pressure corr ocean | `0.137` | `0.208` |
| pressure zonal-anomaly corr all | `0.568` | `0.577` |
| pressure zonal-anomaly corr land | `0.578` | `0.585` |
| pressure zonal-anomaly corr ocean | `0.563` | `0.571` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |

Next M2 step:

1. Keep v4 causal source diagnostics and object extraction.
2. Work on final-pressure smoothness: convert the causal source objects into
   stronger but still bounded stationary-wave pressure geometry.
3. Inspect Earth residuals for where the final field still misses Aleutian /
   Icelandic displacement, Siberian/North American highs, and Southern Ocean
   wave amplitude.
4. Continue Earth-only pressure replay; R2b/R3/R4 remain blocked.

## M2 Pressure-Genesis v6 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v6_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v4:

- Added bounded source-to-pressure wave transfer archived as
  `atmosphere.pressure_genesis_wave_transfer`.
- The transfer projects M2 causal sources into the final pressure proxy while
  keeping R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen,
  biomes, and generated worlds frozen.
- A non-polar land gate prevents Antarctic ice-cap edge structure from being
  amplified as ordinary continental stationary-wave transfer.
- Existing causal support archives remain active:
  `atmosphere.pressure_genesis_source`,
  `atmosphere.ocean_pressure_low_source_support`,
  `atmosphere.ocean_pressure_high_source_support`,
  `atmosphere.land_pressure_source_support`, and
  `atmosphere.terrain_pressure_wave_source_support`.

Map read:

- M0 is no longer the blocker: Atlantic, Pacific, Indian, Arctic, and Southern
  Ocean support remains readable.
- M1 is readable: land/ocean heat-capacity contrast, seasonal land thermal
  anomaly, mixed-layer thermal anomaly, SST-front support, same-latitude SST
  anomaly, lapse cooling, and snow/ice support are all present for attribution.
- DJF North Pacific and North Atlantic causal low sources are visible in the
  source panel.  JJA Southern Ocean causal source is sectorized rather than one
  uniform annular band.
- The v6 wave-transfer panel shows that the previous Antarctic polar-edge
  transfer artifact is reduced.  However, continental transfer still mostly
  follows broad thermal blocks, and the final pressure proxy remains smoother
  than Earth SLP.
- R2a remains unaccepted.  The next missing mechanism is not another global
  amplitude knob; it is source-to-pressure spatial propagation organized by
  terrain, coast, SST-front, and basin waveguides.

Metrics versus the M0 basin checkpoint:

| Metric | M0 basin checkpoint | M2 v6 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.318` |
| standardized-pressure MAE land | `not recorded` | `0.419` |
| standardized-pressure MAE ocean | `0.276` | `0.277` |
| standardized-pressure corr all | `0.539` | `0.548` |
| standardized-pressure corr land | `0.668` | `0.670` |
| standardized-pressure corr ocean | `0.137` | `0.209` |
| pressure zonal-anomaly corr all | `0.568` | `0.579` |
| pressure zonal-anomaly corr land | `0.578` | `0.589` |
| pressure zonal-anomaly corr ocean | `0.563` | `0.570` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `not recorded` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.47s`).

Next M2 step:

1. Keep v6 source and wave-transfer diagnostics.
2. Add a bounded propagation kernel or object graph that lets Aleutian /
   Icelandic / Southern Ocean / continental thermal sources spread into
   pressure centers along plausible terrain/coast/SST-front waveguides instead
   of broad isotropic land blobs.
3. Re-render the same Earth-only pressure packet after each change and read the
   Earth SLP / replay / residual / source / transfer maps before reading
   metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v58 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v58_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/pressure_standardized_delta_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/replay_land_pressure_source_support_seasons.png`

Implementation delta from v55:

- Added DJF Atlantic-Arctic gateway low support and DJF North Pacific /
  Aleutian low support in M2 pressure genesis.  Both are object terms gated by
  semantic basin, latitude/longitude, shelf support, and SST-front support.
- Strengthened the Atlantic-Arctic branch after v56 showed that Greenland Sea
  and Barents/Kara remained under-deep while Aleutian was already close.
- Added DJF North America winter-high relief support.  This is a negative
  land-pressure source over North America, gated by low elevation,
  coast-strength, low interiority, and terrain shelter, to represent maritime
  erosion of the smaller North American winter high compared with Siberia.
- No R2b wind, R3 current, SST, precipitation, biome, or generated-world
  fitting change was made.

Map read:

- DJF Icelandic/Nordic residual improves `+0.188 -> +0.043`.
- DJF Greenland Sea residual improves `+0.400 -> +0.145`; this remains a
  small under-deepness, but the previous winter North Atlantic / Arctic gateway
  blocker is no longer dominant.
- DJF Barents/Kara improves `+0.216 -> +0.051`; Labrador/Baffin improves
  `+0.158 -> +0.038`.
- DJF Aleutian improves `+0.112 -> +0.050`; Bering/Chukchi changes
  `-0.048 -> -0.092`, which is bounded and not a severe over-deep tradeoff.
- DJF Canada land over-high improves `+0.425 -> +0.244`; Alaska improves
  `+0.402 -> +0.301`; Siberia stays near target (`-0.001 -> +0.023`).
- MAM, JJA, and SON are intentionally unchanged by this checkpoint.  Remaining
  owners are MAM Canada land under-high, MAM Arctic cap under-high, JJA Gulf
  Alaska under-high, and broader SON residual structure.

Metrics versus v55:

| Metric | M2 v55 | M2 v58 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.278` | `0.277` |
| standardized-pressure MAE land | `0.387` | `0.385` |
| standardized-pressure MAE ocean | `0.234` | `0.234` |
| standardized-pressure corr all | `0.635` | `0.638` |
| standardized-pressure corr land | `0.688` | `0.690` |
| standardized-pressure corr ocean | `0.529` | `0.534` |
| pressure zonal-anomaly corr all | `0.607` | `0.610` |
| pressure zonal-anomaly corr land | `0.600` | `0.603` |
| pressure zonal-anomaly corr ocean | `0.623` | `0.626` |
| DJF standardized-pressure MAE | `0.227` | `0.222` |
| DJF zonal-anomaly corr | `0.713` | `0.720` |
| MAM standardized-pressure MAE | `0.296` | `0.296` |
| JJA standardized-pressure MAE | `0.203` | `0.203` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.94s`).

Next M2 step:

1. Keep v58 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: MAM Canada / Arctic cap
   under-high, JJA Gulf Alaska under-high, or SON pressure-wave residuals.
3. Preserve v49-v55 fixes plus v58 DJF Atlantic-Arctic / Aleutian lows and
   North America winter-high relief.
4. R2b wind, R3 currents, SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v59 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v59_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/pressure_standardized_delta_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/replay_ocean_pressure_high_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/replay_land_pressure_source_support_seasons.png`

Implementation delta from v58:

- Added a MAM central/eastern North America plains high object.  The support
  is deliberately restricted away from Alaska/Yukon and uses latitude/
  longitude, low-elevation plains/shield terrain, low interiority,
  coast-strength memory, and terrain shelter.
- Added MAM Canadian Arctic and central-west Arctic freeze-high objects.  These
  use semantic Arctic/Atlantic basin support, near-freezing SST, shelf/deep
  basin weighting, SST-front support, and explicit Barents/Kara protection.
- No R2b wind, R3 current, SST, precipitation, biome, or generated-world
  fitting change was made.

Map read:

- MAM Canada land residual improves `-0.309 -> -0.015`.
- MAM west Canada improves `-0.223 -> -0.012`; central Canada improves
  `-0.570 -> -0.261`; east Canada improves `-0.275 -> +0.009`; lower North
  America improves `-0.298 -> +0.030`.
- MAM Alaska/Yukon is not worsened (`+0.258 -> +0.226`), so the new plains
  object avoids the prior broad North America tradeoff.
- MAM Arctic cap improves `-0.289 -> -0.231`; Beaufort improves
  `-0.266 -> -0.068`; Canadian Archipelago improves `-0.475 -> -0.255`.
- MAM Barents/Kara improves `+0.097 -> +0.058`; the new Arctic support does
  not introduce a new Barents over-high wall.
- Baffin/Labrador remains under-high (`-0.198 -> -0.199`) and central Arctic
  remains under-high (`-0.345 -> -0.278`).  These are possible next owners.

Metrics versus v58:

| Metric | M2 v58 | M2 v59 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.277` | `0.276` |
| standardized-pressure MAE land | `0.385` | `0.381` |
| standardized-pressure MAE ocean | `0.234` | `0.233` |
| standardized-pressure corr all | `0.638` | `0.640` |
| standardized-pressure corr land | `0.690` | `0.691` |
| standardized-pressure corr ocean | `0.534` | `0.539` |
| pressure zonal-anomaly corr all | `0.610` | `0.613` |
| pressure zonal-anomaly corr land | `0.603` | `0.607` |
| pressure zonal-anomaly corr ocean | `0.626` | `0.626` |
| MAM standardized-pressure MAE | `0.296` | `0.290` |
| MAM zonal-anomaly corr | `0.604` | `0.615` |
| DJF standardized-pressure MAE | `0.222` | `0.222` |
| JJA standardized-pressure MAE | `0.203` | `0.203` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.31s`).

Next M2 step:

1. Keep v59 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: residual MAM Baffin-Labrador /
   central Arctic under-high, JJA Gulf Alaska / Aleutian under-high, or SON
   pressure-wave residuals.
3. Preserve v49-v59 pressure-source fixes and keep downstream systems blocked.

## M2 Pressure-Genesis v60 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v60_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/pressure_standardized_delta_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/replay_ocean_pressure_high_source_support_seasons.png`

Implementation delta from v59:

- Added a JJA eastern / central North Pacific high-pressure object in M2
  pressure genesis.
- The support is restricted to semantic Pacific basin cells in the Gulf Alaska
  / eastern Aleutian / central North Pacific sector and is gated by latitude,
  longitude, shelf/front support, cool same-latitude SST, and an Arctic taper.
- The support suppresses the already over-high NW Pacific 130E-160E sector and
  does not touch Beaufort/Arctic or North Atlantic pressure support.
- No R2b wind, R3 current, SST, precipitation, biome, or generated-world
  fitting change was made.

Map read:

- JJA Gulf Alaska residual improves `-0.265 -> -0.118`.
- JJA Gulf Alaska east / west improve `-0.235 -> -0.089` and
  `-0.309 -> -0.170`.
- JJA eastern Aleutian improves `-0.343 -> -0.181`; western Aleutian changes
  `-0.101 -> +0.014`, which remains bounded.
- JJA North Pacific 160W-130W and 180-160W improve `-0.247 -> -0.118` and
  `-0.251 -> -0.136`.
- NW Pacific 130E-160E remains over-high but does not worsen
  (`+0.369 -> +0.360`); Beaufort remains protected (`+0.303 -> +0.299`).
- JJA North Atlantic subpolar remains under-high (`-0.216 -> -0.228`) and is
  a possible later owner.

Metrics versus v59:

| Metric | M2 v59 | M2 v60 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.276` | `0.275` |
| standardized-pressure MAE land | `0.381` | `0.381` |
| standardized-pressure MAE ocean | `0.233` | `0.232` |
| standardized-pressure corr all | `0.640` | `0.642` |
| standardized-pressure corr ocean | `0.539` | `0.543` |
| pressure zonal-anomaly corr all | `0.613` | `0.615` |
| pressure zonal-anomaly corr ocean | `0.626` | `0.631` |
| JJA standardized-pressure MAE | `0.203` | `0.199` |
| JJA zonal-anomaly corr | `0.741` | `0.749` |
| MAM standardized-pressure MAE | `0.290` | `0.290` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.41s`).

Next M2 step:

1. Keep v60 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: residual MAM Baffin-Labrador /
   central Arctic under-high, JJA North Atlantic subpolar under-high, or SON
   pressure-wave residuals.
3. Preserve v49-v60 pressure-source fixes and keep downstream systems blocked.

## M2 Pressure-Genesis v9 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v9_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v8:

- Added directional masked diffusion helpers for M2 pressure genesis.
- Ocean source-to-pressure diffusion now uses a storm-track/SST-front axis
  instead of only scalar conductance.
- Land source-to-pressure diffusion now uses terrain-barrier and coastal axes
  instead of only scalar conductance.
- The v8 object projection and non-polar land gate remain in place.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- v9 preserves v8's DJF North Pacific / North Atlantic low-pressure transfer
  footprints and Eurasian winter-high footprint.
- JJA Southern Ocean remains segmented and wave-like.  No new Antarctic
  polar-edge artifact is visible.
- The final pressure proxy still looks too similar to v8 and remains smoother
  than Earth SLP.  The transfer/source diagnostics are more physically
  organized than the final pressure response.
- R2a remains unaccepted.  The remaining M2 owner is final-pressure expression
  and placement strength: pressure-center objects need to affect the pressure
  proxy more coherently without global amplitude painting.

Metrics versus v8:

| Metric | M2 v8 | M2 v9 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE land | `0.418` | `0.418` |
| standardized-pressure MAE ocean | `0.277` | `0.277` |
| standardized-pressure corr all | `0.550` | `0.550` |
| standardized-pressure corr land | `0.672` | `0.672` |
| standardized-pressure corr ocean | `0.212` | `0.212` |
| pressure zonal-anomaly corr all | `0.581` | `0.581` |
| pressure zonal-anomaly corr land | `0.592` | `0.592` |
| pressure zonal-anomaly corr ocean | `0.570` | `0.570` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.57s`).

Next M2 step:

1. Keep v9 directional diffusion and v8 object projection.
2. Stop treating anisotropy alone as the blocker; it is now present but too
   weak to solve the final pressure map.
3. Improve how projected pressure-center objects enter the final standardized
   pressure proxy.  Candidate mechanism: bounded object-center pressure solve
   that preserves the broad thermal base while increasing source-supported
   Aleutian/Icelandic/Southern Ocean/continental center contrast.
4. Continue Earth-only pressure replay; R2b/R3/R4/R5/R6/R7/R8 and generated
   worlds remain blocked as fitting targets.

## M2 Pressure-Genesis v10 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v10_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v9:

- Added a bounded final-pressure expression solve inside M2.
- The expression solve starts from the provisional final pressure after source
  and directional transfer, extracts nonzonal pressure anomaly, and enhances it
  only where source support is present and the anomaly sign aligns with the
  signed M2 source/transfer field.
- The expression increment is latitude-band centered and clipped, so it should
  not become a global high-latitude band.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- DJF transfer now expresses North Pacific / North Atlantic lows and the
  Eurasian winter high more forcefully than v9.
- JJA Southern Ocean still reads as a segmented wave structure, not a single
  annular paint band.
- The final standardized pressure proxy is improved but still too smooth
  versus Earth SLP.  Pressure-center structure is now more visible in final
  pressure than in v9, but not yet realistic enough for R2a promotion.
- R2a remains unaccepted.  Remaining M2 owner: center placement and shape in
  final pressure, especially making Aleutian/Icelandic/Southern Ocean centers
  coherent without over-amplifying high-latitude bands.

Metrics versus v9:

| Metric | M2 v9 | M2 v10 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE land | `0.418` | `0.416` |
| standardized-pressure MAE ocean | `0.277` | `0.276` |
| standardized-pressure corr all | `0.550` | `0.550` |
| standardized-pressure corr land | `0.672` | `0.671` |
| standardized-pressure corr ocean | `0.212` | `0.220` |
| pressure zonal-anomaly corr all | `0.581` | `0.582` |
| pressure zonal-anomaly corr land | `0.592` | `0.593` |
| pressure zonal-anomaly corr ocean | `0.570` | `0.571` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.19s`).

Next M2 step:

1. Keep v10 final-pressure expression as the current baseline.
2. Improve placement/shape rather than amplitude: centers need better spatial
   coherence in final pressure, especially Aleutian/Icelandic lows and
   Southern Ocean waves.
3. Add a guard or diagnostic for high-latitude expression becoming too
   band-like before increasing expression strength.
4. Continue Earth-only pressure replay; R2b/R3/R4/R5/R6/R7/R8 and generated
   worlds remain blocked as fitting targets.

## M2 Pressure-Genesis v8 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v8_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v7:

- Added bounded object-level pressure-center projection inside M2.
- The projection uses already-audited causal supports:
  ocean-low support, land thermal-center source support, and terrain-wave
  source support.
- The v7 waveguide transfer remains the baseline but is stronger, so projected
  objects can affect the final pressure proxy rather than staying nearly
  invisible in diagnostics.
- The non-polar land gate remains in place; Antarctic ice-cap edge structure is
  not treated as ordinary continental stationary-wave forcing.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- DJF transfer now shows clearer North Pacific / North Atlantic low-pressure
  footprints and a stronger Eurasian winter high footprint than v7.
- JJA Southern Ocean remains segmented and wave-like rather than a single
  uniform annular band.  No obvious Antarctic polar-edge artifact returned.
- The final pressure proxy is still visibly smoother than Earth SLP.  The
  pressure centers are more readable in transfer than in the final standardized
  pressure map, so R2a is still not accepted.
- Remaining owner: M2 pressure-center placement and footprint geometry.  The
  next step should improve anisotropic source projection along storm-track,
  terrain, coast, and SST-front supports instead of increasing global
  amplitude.

Metrics versus v7:

| Metric | M2 v7 | M2 v8 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE land | `0.418` | `0.418` |
| standardized-pressure MAE ocean | `0.277` | `0.277` |
| standardized-pressure corr all | `0.549` | `0.550` |
| standardized-pressure corr land | `0.671` | `0.672` |
| standardized-pressure corr ocean | `0.209` | `0.212` |
| pressure zonal-anomaly corr all | `0.581` | `0.581` |
| pressure zonal-anomaly corr land | `0.591` | `0.592` |
| pressure zonal-anomaly corr ocean | `0.571` | `0.570` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 70.93s`).

Next M2 step:

1. Keep v8 object projection and v7 waveguide diffusion as the baseline.
2. Improve pressure-center placement and footprint anisotropy: use source
   object centroids, coast/terrain orientation, SST-front orientation, and
   basin/continent labels to elongate Aleutian/Icelandic/Southern Ocean/
   continental centers in plausible directions instead of isotropic blobs.
3. Re-render the same Earth-only pressure packet and read Earth SLP, replay
   pressure, residual, source, transfer, and M0/M1 support maps before reading
   metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v7 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v7_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v6:

- Added a local weighted-diffusion helper for masked fields.
- Replaced the v6 mostly isotropic source-to-pressure transfer with a
  geography-weighted waveguide transfer.
- Ocean transfer is weighted by SST-front support, same-latitude SST anomaly,
  open-ocean exposure, and subpolar/basin exposure.
- Land transfer is weighted by coast-strength, terrain barrier support,
  land-source support, and terrain-wave support, with the same non-polar gate
  that prevents Antarctic ice-cap edge amplification.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- The M2 source map is intentionally unchanged from v6: DJF North Pacific /
  North Atlantic low sources and JJA Southern Ocean source sectors remain
  visible and auditable.
- The v7 transfer map is more structured along the Northern Hemisphere
  storm-track/coastal waveguide than v6 and does not reintroduce the Antarctic
  polar-edge artifact.
- The final pressure replay is still close to v6: it remains too smooth and
  too block-like relative to Earth SLP.  The replay still lacks sufficiently
  strong, readable Aleutian / Icelandic / Southern Ocean pressure-center
  geometry and still underuses terrain/coast/SST-front waveguides in the final
  pressure map.
- R2a remains unaccepted.  v7 is a correct mechanism-direction step, but the
  next owner is stronger object-level source-to-pressure projection, not R2b
  wind or downstream climate tuning.

Metrics versus v6:

| Metric | M2 v6 | M2 v7 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.318` | `0.317` |
| standardized-pressure MAE land | `0.419` | `0.418` |
| standardized-pressure MAE ocean | `0.277` | `0.277` |
| standardized-pressure corr all | `0.548` | `0.549` |
| standardized-pressure corr land | `0.670` | `0.671` |
| standardized-pressure corr ocean | `0.209` | `0.209` |
| pressure zonal-anomaly corr all | `0.579` | `0.581` |
| pressure zonal-anomaly corr land | `0.589` | `0.591` |
| pressure zonal-anomaly corr ocean | `0.570` | `0.571` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.79s`).

Next M2 step:

1. Keep v7 waveguide diffusion as the local propagation baseline.
2. Add object-level pressure-center projection: use the existing low/high/
   land/terrain source supports to create bounded Aleutian, Icelandic,
   Southern Ocean, Siberian/North American winter-high, and continental
   summer-low pressure-center objects with explicit footprint constraints.
3. The acceptance packet remains the same Earth-only pressure contact sheet.
   Read Earth SLP, replay pressure, residual, source, transfer, and M0/M1
   support maps before reading metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v16 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v16_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v10:

- v11/v12 object-local final-pressure expression was rejected: it made the map
  no more readable and degraded pressure metrics.
- v13/v14 stronger object/source expression was rejected for the same reason.
- v16 keeps the v10 final-pressure expression solve and adds one source-side
  mechanism.  Winter cold continental highs are spread downwind along the
  westerly/storm-track axis into nearby open ocean, then used as support for
  subpolar ocean-low source selection.  Southern Hemisphere downwind support is
  damped so the Southern Ocean front/wavenumber gate remains the primary
  Southern Ocean pressure-wave source.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- DJF North Pacific / North Atlantic source placement remains readable and is
  now partly derived from cold-continent downstream geometry rather than only
  broad open-ocean latitude/front scoring.
- JJA Southern Ocean source support remains segmented; the new mechanism did
  not reintroduce an Antarctic edge artifact or a full circumpolar ring.
- The final standardized pressure map is still too smooth relative to Earth
  SLP.  v16 is retained as a small mechanism-direction step, not an R2a
  acceptance point.

Metrics versus v10:

| Metric | M2 v10 | M2 v16 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.317` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.276` | `0.276` |
| standardized-pressure corr all | `0.550` | `0.551` |
| standardized-pressure corr land | `0.671` | `0.671` |
| standardized-pressure corr ocean | `0.220` | `0.221` |
| pressure zonal-anomaly corr all | `0.582` | `0.583` |
| pressure zonal-anomaly corr land | `0.593` | `0.594` |
| pressure zonal-anomaly corr ocean | `0.571` | `0.572` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.24s`).

Next M2 step:

1. Keep v16 downwind cold-continent support only as source-side assistance; do
   not increase it into a hard-coded pressure painting term.
2. The remaining blocker is final-pressure dominance by the smoothed upstream
   pressure proxy.  Re-balance the source-supported M2 pressure field against
   that upstream proxy while preserving v10's bounded nonzonal expression.
3. Continue Earth-only map-read acceptance: Earth SLP, replay pressure,
   residual, M2 source, M2 transfer, M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v17 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v17_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v16:

- Reconstructed the v16 pressure field as `upstream_base + M2_source +
  M2_transfer` before editing code.
- Probe result: direct raw-source strengthening makes the SLP pattern worse,
  while stronger source-to-pressure transfer improves the Earth replay.
- v17 keeps v16 source geometry, uses the full source field as a causal driver,
  reduces direct source pressure contribution to `0.80`, and increases the
  bounded transfer contribution to `1.45`.  The archived
  `atmosphere.pressure_genesis_wave_transfer` now represents the stronger
  transfer contribution used by final pressure.
- R2b wind, R3 currents, SST, moisture, precipitation, ice, Koppen, biomes, and
  generated worlds remain unchanged as fitting targets.

Map read:

- The M2 transfer panel is visibly stronger in Northern Hemisphere subpolar
  source-supported regions and in the JJA Southern Ocean wave belt.
- Final pressure is incrementally less dominated by the smoothed upstream
  proxy, especially in JJA.
- Final pressure centers are still too smooth and insufficiently coherent
  relative to Earth SLP.  v17 is retained as a synthesis-layer improvement, not
  an R2a acceptance point.

Metrics versus v16:

| Metric | M2 v16 | M2 v17 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.317` | `0.316` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.276` | `0.275` |
| standardized-pressure corr all | `0.551` | `0.552` |
| standardized-pressure corr land | `0.671` | `0.672` |
| standardized-pressure corr ocean | `0.221` | `0.222` |
| pressure zonal-anomaly corr all | `0.583` | `0.584` |
| pressure zonal-anomaly corr land | `0.594` | `0.595` |
| pressure zonal-anomaly corr ocean | `0.572` | `0.573` |
| DJF pressure zonal-anomaly corr | `0.710` | `0.709` |
| JJA pressure zonal-anomaly corr | `0.723` | `0.725` |
| SON pressure zonal-anomaly corr | `0.373` | `0.371` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 104.60s`).

Next M2 step:

1. Keep the v17 source/transfer balance unless visual map reading shows a
   stronger reason to revert it.
2. Improve pressure-center morphology and placement inside M2: centers should
   become more coherent in final pressure without hard-coded seasonal gain
   tables or broad latitude-band painting.
3. Continue Earth-only map-read acceptance: Earth SLP, replay pressure,
   residual, M2 source, M2 transfer, M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v18 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v18_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v17:

- Kept v17's source/transfer balance: source remains a causal trigger, while
  transfer carries more of the actual pressure expression.
- Added a bounded morphology adjustment to the already-projected transfer:
  ocean-low support deepens negative transfer cores; ocean-high, land-source,
  and terrain-wave supports damp broad or fragmented transfer expression.
- This is still M2-only.  It uses existing support fields and does not add
  hard-coded coordinates, seasonal gain tables, downstream wind/current/SST
  tuning, precipitation tuning, or biome fitting.

Map read:

- The transfer map is slightly more center-like and less terrain-fragmented than
  v17.
- No new Antarctic edge artifact, Southern Ocean ring, or high-latitude speckle
  appears in the pressure replay.
- Final pressure is still too smooth and weakly center-organized compared with
  Earth SLP.  v18 is retained as a balanced morphology checkpoint, not an R2a
  acceptance point.

Metrics versus v17:

| Metric | M2 v17 | M2 v18 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.316` | `0.316` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.275` | `0.275` |
| standardized-pressure corr all | `0.552` | `0.553` |
| standardized-pressure corr land | `0.672` | `0.673` |
| standardized-pressure corr ocean | `0.222` | `0.224` |
| pressure zonal-anomaly corr all | `0.584` | `0.585` |
| pressure zonal-anomaly corr land | `0.595` | `0.596` |
| pressure zonal-anomaly corr ocean | `0.573` | `0.574` |
| DJF pressure zonal-anomaly corr | `0.709` | `0.709` |
| JJA pressure zonal-anomaly corr | `0.725` | `0.724` |
| SON pressure zonal-anomaly corr | `0.371` | `0.371` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 204.10s`).

Next M2 step:

1. Keep v18 as the current balanced checkpoint unless visual map reading shows a
   regression.
2. The remaining R2a blocker is final-pressure center coherence and placement,
   not source-map readability or downstream wind/current behavior.
3. Continue Earth-only map-read acceptance: Earth SLP, replay pressure,
   residual, M2 source, M2 transfer, M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v19 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v19_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v18:

- Kept v18's source/transfer balance and bounded transfer-morphology adjustment.
- Added a bounded thermal-phase adjustment to the already-projected transfer.
  Same-latitude SST and land-temperature anomalies now decide whether
  high-support, low-support, and land-source footprints strengthen or weaken
  final pressure.  This addresses the v18 attribution where high support had
  the wrong seasonal sign over warm subtropical ocean in MAM while cold phases
  needed stronger high expression in SON.
- This is still M2-only.  It uses existing M1/M2 fields and does not add
  coordinate-specific targets, hard-coded seasonal gain tables, downstream
  wind/current/SST tuning, precipitation tuning, or biome fitting.

Map read:

- The pressure-transfer panel is more seasonally differentiated than v18,
  especially across MAM/SON.  Final pressure also picks up a little more
  ocean-center phase structure.
- No full Southern Ocean ring, Antarctic edge artifact, or heat-wall latitude
  band is visible in the contact sheet.
- Final pressure remains too smooth and not center-organized enough compared
  with Earth SLP.  v19 is retained as the current checkpoint, not an R2a
  acceptance point.
- Watch item: thermal-phase support creates broader subtropical transfer
  patches.  The next M2 pass should sharpen pressure-center morphology and
  placement rather than simply increasing amplitude.

Metrics versus v18:

| Metric | M2 v18 | M2 v19 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.316` | `0.310` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.275` | `0.268` |
| standardized-pressure corr all | `0.553` | `0.564` |
| standardized-pressure corr land | `0.673` | `0.675` |
| standardized-pressure corr ocean | `0.224` | `0.269` |
| pressure zonal-anomaly corr all | `0.585` | `0.587` |
| pressure zonal-anomaly corr land | `0.596` | `0.600` |
| pressure zonal-anomaly corr ocean | `0.574` | `0.574` |
| DJF pressure zonal-anomaly corr | `0.709` | `0.708` |
| MAM pressure zonal-anomaly corr | `0.556` | `0.582` |
| JJA pressure zonal-anomaly corr | `0.724` | `0.723` |
| SON pressure zonal-anomaly corr | `0.371` | `0.371` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 139.05s`).

Next M2 step:

1. Keep v19 as the current checkpoint, but do not promote R2a.
2. The remaining blocker is final-pressure center coherence and placement:
   Aleutian/Icelandic/Southern Ocean and continental seasonal centers must
   become coherent in final pressure, not only in support/source diagnostics.
3. Avoid amplitude-only tuning.  The next mechanism should sharpen object
   footprints and center placement using M0 basin/coast/barrier supports and M1
   energy/frontal phase.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v25 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v25_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v19:

- Tightened the Southern Ocean pressure-wave gate by reducing its sector floor
  and sharpening the wave core.  This reduces continuous ring-like source
  expression while preserving the wavenumber/front basis.
- Changed the final pressure expression from a uniform `0.80 * source` term to
  a domain-weighted direct-source term.  Land and Southern Ocean source
  expression are reduced; North Hemisphere mid/high-latitude ocean-low cores
  retain more direct expression so Aleutian/Icelandic lows are not erased.
- Rejected trials: v20/v21 allowed shoulder-season low support too broadly and
  produced unrealistic MAM/SON ocean-low fields; v23 basin-area reweighting
  did not improve map read and was removed.
- This is still M2-only.  It uses existing M0/M1/M2 supports and does not add
  coordinate-specific targets, hard-coded seasonal gain tables, downstream
  wind/current/SST tuning, precipitation tuning, or biome fitting.

Map read:

- JJA Southern Ocean over-deep residual is reduced relative to v19, and the
  source/transfer panels remain free of a full Southern Ocean ring or Antarctic
  edge artifact.
- DJF Aleutian Low remains close to v19.  The uniform source reduction trial
  weakened it too much; v25's domain weighting avoids that regression.
- DJF Icelandic Low remains too weak.  SON North Pacific still has almost no
  low-source support.  Final pressure is therefore still too smooth and weakly
  center-organized compared with Earth SLP.
- v25 is retained as the current checkpoint, not an R2a acceptance point.

Metrics versus v19:

| Metric | M2 v19 | M2 v25 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.310` | `0.309` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.268` | `0.266` |
| standardized-pressure corr all | `0.564` | `0.564` |
| standardized-pressure corr land | `0.675` | `0.674` |
| standardized-pressure corr ocean | `0.269` | `0.272` |
| pressure zonal-anomaly corr all | `0.587` | `0.589` |
| pressure zonal-anomaly corr land | `0.600` | `0.600` |
| pressure zonal-anomaly corr ocean | `0.574` | `0.577` |
| DJF pressure zonal-anomaly corr | `0.708` | `0.707` |
| MAM pressure zonal-anomaly corr | `0.582` | `0.585` |
| JJA pressure zonal-anomaly corr | `0.723` | `0.726` |
| SON pressure zonal-anomaly corr | `0.371` | `0.374` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 115.70s`).

Next M2 step:

1. Keep v25 as the current checkpoint, but do not promote R2a.
2. The remaining blocker is source placement, not downstream response:
   strengthen/relocate the DJF Icelandic Low using basin/coast/front/terrain
   supports, and create a physically justified SON North Pacific source without
   reactivating broad MAM/SON ocean-low bands.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v30 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v30_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v25:

- Added a post-normalized Northern Hemisphere winter coastal/front supplement.
  This targets the v25 Icelandic Low attribution: strong SST-front support
  existed, but the source score over-suppressed basin-edge frontal lows through
  the `open_ocean` gate.
- Added a bounded shoulder-season warm-ocean/front low-source candidate.  The
  candidate uses positive seasonal SST anomaly, SST-front support, open-ocean
  exposure, subpolar latitude, and basin component selection.  After object
  selection, its amplitude is capped and modulated by open-ocean exposure so it
  cannot normalize to a winter-strength low.
- Rejected trials: v26 improved Icelandic but weakened Aleutian by changing the
  global low-support normalization; v28 over-deepened SON ocean lows; v29 was
  better but still too strong.  v30 keeps the mechanism with lower shoulder
  amplitude.
- This is still M2-only.  It uses existing M0/M1/M2 supports and does not add
  coordinate-specific targets, hard-coded seasonal gain tables, downstream
  wind/current/SST tuning, precipitation tuning, or biome fitting.

Map read:

- SON North Pacific now has a weak local low-source/transfer expression.  v25
  had essentially none.
- MAM North Pacific/North Atlantic remain untriggered, so the shoulder-season
  rule is not just a symmetric equinox low-source paint.
- DJF Icelandic Low is slightly stronger than v25; DJF Aleutian remains close
  to v25.
- Final pressure remains too smooth.  The Icelandic center is still weak, and
  the new SON source needs better final-pressure expression without damaging
  ocean/zonal correlations.  v30 is retained as the current checkpoint, not an
  R2a acceptance point.

Metrics versus v25:

| Metric | M2 v25 | M2 v30 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.309` | `0.309` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.266` | `0.266` |
| standardized-pressure corr all | `0.564` | `0.563` |
| standardized-pressure corr land | `0.674` | `0.674` |
| standardized-pressure corr ocean | `0.272` | `0.268` |
| pressure zonal-anomaly corr all | `0.589` | `0.585` |
| pressure zonal-anomaly corr land | `0.600` | `0.598` |
| pressure zonal-anomaly corr ocean | `0.577` | `0.569` |
| DJF pressure zonal-anomaly corr | `0.707` | `0.707` |
| MAM pressure zonal-anomaly corr | `0.585` | `0.585` |
| JJA pressure zonal-anomaly corr | `0.726` | `0.726` |
| SON pressure zonal-anomaly corr | `0.374` | `0.370` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.16s`).

Next M2 step:

1. Keep v30 as a historical source-placement checkpoint, but do not promote R2a.
2. The remaining blocker is final-pressure expression and compactness: preserve
   the new SON North Pacific source while recovering ocean/zonal correlations,
   and strengthen the Icelandic Low without weakening Aleutian.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v38 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v38_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v30:

- The shoulder-season warm-ocean/front source remains M2-only, but its
  post-object amplitude now carries stronger basin-scale and open-ocean
  weighting.  This prevents a small/near-coastal basin front from normalizing
  to the same final low strength as a broad open-ocean storm-track sector.
- A winter-only Northern Hemisphere subpolar source-expression boost deepens
  existing Aleutian/Icelandic-type lows without changing MAM/SON/JJA source
  activation.
- Southern Ocean source placement no longer uses the fixed longitude wave as
  the dominant gate.  Front support, shelf/slope support, and same-latitude SST
  anomaly dominate; the longitude wave is only a weak perturbation.
- A signed Southern Ocean wave-transfer term converts high-support sectors into
  relative lows and low-support sectors into relative highs.  This was added
  because v34/v36 source maps could change low-source strength, but could not
  create the positive/negative pressure wave visible in Earth SLP.
- Rejected intermediate trials: v31 erased the useful SON North Pacific source;
  v32/v33 restored it and reduced the North Atlantic excess; v35/v36 removed
  too much Southern Ocean source; v37 introduced the correct signed transfer
  direction; v38 increases the wave amplitude enough to be visible.
- This is still M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- DJF Aleutian and Icelandic lows are modestly deeper than v30.  Region means
  improve from v33 to v34/v38 without weakening Aleutian.
- SON North Pacific keeps a weak local low source.  SON North Atlantic remains
  active, but the v30 over-deepening is reduced.
- JJA Southern Ocean now has a signed transfer wave.  In the coarse
  `-65..-45` latitude band read, the `60..120E` sector moves from replay
  `+0.084` in v34 to `-0.006` in v38 against Earth `-0.173`; false-low
  sectors such as `-180..-120` are lifted from `-0.019` to `+0.088` against
  Earth `+0.130`.
- Remaining defect: the final pressure field is still too smooth and not yet
  center-organized enough.  Southern Ocean phase is improved but still weak,
  and North Atlantic/Icelandic placement remains incomplete.  R2a is not
  accepted.

Metrics versus v30:

| Metric | M2 v30 | M2 v38 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.309` | `0.309` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.266` | `0.266` |
| standardized-pressure corr all | `0.563` | `0.564` |
| standardized-pressure corr land | `0.674` | `0.674` |
| standardized-pressure corr ocean | `0.268` | `0.272` |
| pressure zonal-anomaly corr all | `0.585` | `0.588` |
| pressure zonal-anomaly corr land | `0.598` | `0.600` |
| pressure zonal-anomaly corr ocean | `0.569` | `0.576` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.11s`).
- Full test suite was started and interrupted for time after
  `103 passed in 444.62s`; no failure was observed before interruption.

Next M2 step:

1. Keep v38 as a historical Southern Ocean wave checkpoint, but do not promote R2a.
2. The remaining blocker is pressure-center compactness and placement:
   continue tightening Icelandic/North Atlantic and Southern Ocean pressure
   objects without using downstream wind/current/SST/precipitation/biome
   fitting.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v40 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v40_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/replay_m0_boundary_support_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v38:

- Added a bounded Northern Hemisphere winter subpolar low-core expression term
  in the M2 source-to-pressure transfer.  It uses existing `low_support`,
  SST-front support, winter hemisphere gating, and the North-Hemisphere
  subpolar ocean mask.
- This is not a new source generator.  It only lets already-audited
  Aleutian/Icelandic-type source cores express more compactly in final
  pressure.
- v39 tested the same direction with lower amplitude.  Offline replay on v39
  arrays showed a coefficient near `0.20` improved Aleutian/Icelandic residuals
  and ocean pressure correlation before global MAE began to flatten, so v40
  uses the stronger bounded value.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- DJF Aleutian core replay improves from `-0.575` in v38 to `-0.645` in v40
  against Earth `-0.840`.
- DJF Icelandic core improves from `-0.435` in v38 to `-0.493` in v40 against
  Earth `-0.706`.
- DJF North Pacific east, Labrador edge, and NE Atlantic residuals also move in
  the right direction, but remain under-deepened.
- MAM North Pacific/North Atlantic remain untriggered.  SON North Pacific and
  North Atlantic are unchanged by this winter-only term.  Southern Ocean v38
  signed wave remains present.
- R2a is still not accepted.  The final pressure map is still too smooth and
  not compact/organized enough around North Atlantic edge/Labrador and
  Southern Ocean wave sectors.

Metrics versus v38:

| Metric | M2 v38 | M2 v40 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.309` | `0.309` |
| standardized-pressure MAE land | `0.416` | `0.416` |
| standardized-pressure MAE ocean | `0.266` | `0.265` |
| standardized-pressure corr all | `0.564` | `0.565` |
| standardized-pressure corr ocean | `0.272` | `0.277` |
| pressure zonal-anomaly corr all | `0.588` | `0.589` |
| pressure zonal-anomaly corr ocean | `0.576` | `0.577` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.26s`).

Next M2 step:

1. Keep v40 as a historical checkpoint, but do not promote R2a.
2. The remaining blocker is pressure-center compactness and placement:
   continue tightening North Atlantic edge/Labrador, remaining DJF low-core
   strength, and Southern Ocean wave sectors.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v44 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v44_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v40:

- Added a Northern Hemisphere DJF coastal-land inheritance term.  It spreads
  adjacent ocean-low support into nearby subpolar coastal/island land through
  the existing coast and interiority gates, then adds it as bounded
  source-to-pressure transfer.
- Added a bounded poleward subpolar SST-front / storm-track low-support floor
  between roughly `50-67 N`.  This prevents Labrador/Icelandic lows from being
  erased by the earlier open-ocean exposure gate while keeping the support tied
  to SST-front, basin-scale, winter-hemisphere, and near-coast geometry.
- v41 was too weak; v42 made the coastal inheritance visible.  v43 was rejected
  because putting lee-low support into pressure-object selection overdeepened
  Labrador and NW Pacific while weakening Icelandic low support through object
  competition.  v44 keeps the noncompetitive floor direction.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- DJF Labrador ocean residual improves from `+0.253` in v42 to `+0.179` in
  v44.
- DJF Icelandic ocean residual improves from `+0.353` in v42 to `+0.333` in
  v44.
- DJF Labrador land residual improves from `+0.242` in v42 to `+0.235` in v44.
- DJF Icelandic land residual improves from `+0.512` in v42 to `+0.504` in v44.
- DJF NE Atlantic ocean becomes slightly overdeepened (`-0.014 -> -0.033`), so
  the next pass must avoid broad North Atlantic amplification.
- R2a is still not accepted.  The North Atlantic / Icelandic low remains too
  smooth and weak, Arctic/Nordic winter low structure is still under-expressed,
  and Southern Ocean wave sectors remain incomplete.

Metrics versus v40:

| Metric | M2 v40 | M2 v44 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.309` | `0.308` |
| standardized-pressure MAE land | `0.416` | `0.415` |
| standardized-pressure MAE ocean | `0.265` | `0.265` |
| standardized-pressure corr all | `0.565` | `0.566` |
| standardized-pressure corr ocean | `0.277` | `0.281` |
| pressure zonal-anomaly corr all | `0.589` | `0.590` |
| pressure zonal-anomaly corr ocean | `0.577` | `0.577` |
| DJF standardized-pressure MAE | `0.232` | `0.230` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.40s`).

Next M2 step:

1. Keep v44 as a historical checkpoint, but do not promote R2a.
2. The remaining blocker is residual North Atlantic/Icelandic compactness and
   Southern Ocean wave-sector organization, not downstream wind/current/SST or
   precipitation tuning.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v45 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v45_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v44:

- Added an Atlantic-Arctic gateway low-support floor.  It computes graph
  distance from existing Atlantic subpolar low support through the
  Atlantic/Arctic ocean domain, then applies a bounded floor only to nearby
  `62-80 N` Arctic-basin marginal seas with shelf and SST-front support.
- The term is designed to express North Atlantic pressure influence into the
  Nordic/Greenland/Barents gateway without painting the Pacific-side Arctic.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- DJF Nordic/Arctic gateway residual improves from `+0.642` in v44 to
  `+0.538` in v45.
- DJF Greenland Sea residual improves from `+0.645` to `+0.562`.
- DJF Barents/Kara residual improves from `+0.479` to `+0.425`.
- DJF Icelandic ocean residual improves from `+0.337` to `+0.314`.
- Beaufort/Arctic and Bering/Chukchi remain effectively protected, so v45 does
  not create a spurious pan-Arctic low.
- R2a is still not accepted.  Nordic/Arctic and Icelandic lows remain
  underdeepened, and Southern Ocean shoulder-season / wave-sector organization
  remains the clearest next owner.

Metrics versus v44:

| Metric | M2 v44 | M2 v45 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.308` | `0.308` |
| standardized-pressure MAE land | `0.415` | `0.415` |
| standardized-pressure MAE ocean | `0.265` | `0.265` |
| standardized-pressure corr all | `0.566` | `0.566` |
| standardized-pressure corr ocean | `0.281` | `0.284` |
| pressure zonal-anomaly corr all | `0.590` | `0.590` |
| pressure zonal-anomaly corr ocean | `0.577` | `0.579` |
| DJF standardized-pressure MAE | `0.230` | `0.229` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.25s`).

Next M2 step:

1. Keep v45 as a historical checkpoint, but do not promote R2a.
2. The next owner is Southern Ocean shoulder-season / wave-sector source
   organization.  Residual North Atlantic compactness remains a watch item.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v46 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v46_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v45:

- Added a Southern Ocean shoulder-season low-support floor.  It uses the
  existing Southern Ocean front/shelf/wave gate to create MAM/SON low-source
  support in the semantic Southern Ocean basin.
- This fixes the v45 mechanism gap where MAM/SON Southern Ocean had M1
  front/shelf support but no M2 low-source expression.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- Semantic Southern Ocean MAM residuals improve in every 60-degree sector.
- Semantic Southern Ocean SON residuals improve in all target low-pressure
  sectors except the 0-60E sector, which becomes more overdeepened.
- Wide 45-75S SON residuals improve at `-120..-60`, `60..120`, and
  `120..180`; `-60..60` remains too low.
- R2a is still not accepted.  v46 fixes the missing shoulder-source class, but
  the Southern Ocean response remains too band-like.  The next M2 owner is
  wave-sector sharpening, not broader amplitude.

Metrics versus v45:

| Metric | M2 v45 | M2 v46 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.308` | `0.306` |
| standardized-pressure MAE land | `0.415` | `0.416` |
| standardized-pressure MAE ocean | `0.265` | `0.262` |
| standardized-pressure corr all | `0.566` | `0.571` |
| standardized-pressure corr ocean | `0.284` | `0.308` |
| pressure zonal-anomaly corr all | `0.590` | `0.591` |
| pressure zonal-anomaly corr ocean | `0.579` | `0.577` |
| MAM standardized-pressure MAE | `0.357` | `0.355` |
| SON standardized-pressure MAE | `0.434` | `0.430` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.58s`).

Next M2 step:

1. Keep v46 as the current checkpoint, but do not promote R2a.
2. The next owner is Southern Ocean wave-sector sharpening, not broader
   Southern Ocean amplitude.  Residual North Atlantic compactness remains a
   secondary watch item.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v47/v48 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v47_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v47_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v47_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v48_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v46:

- v47 kept the v46 Southern Ocean shoulder-source class and added M2-only
  shoulder-season source-to-pressure transfer geometry.
- The new transfer uses M0/M1 support already present in the replay:
  semantic Southern Ocean mask, latitude gate, SST-front support, shelf
  support, open-ocean exposure, and same-latitude SST anomaly.
- MAM receives a seasonal phase term that lifts the Indian/Australian
  subantarctic sector and deepens the South Atlantic/east-Pacific storm-track
  side.
- SON receives an Atlantic-sector lift plus a Pacific/Amundsen open-ocean low
  anchor tied to open-ocean and same-latitude SST support.
- v47 was directionally useful for SON but pushed MAM too poleward.  v48 keeps
  the SON repair and narrows the MAM phase support toward the more equatorward
  Southern Ocean storm-track belt, avoiding excessive lifting of the polar-side
  low core.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- v48 changes the MAM/SON transfer map from a mostly band-like shoulder-season
  Southern Ocean response into a readable wave-sector pattern.
- Wide 45-75S SON residuals improve from v46 to v48:
  `-180..-120` goes `+0.71 -> +0.47`, `-120..-60` goes
  `+0.10 -> -0.01`, `0..60` goes `-0.24 -> -0.16`,
  `60..120` goes `+0.29 -> +0.23`, and `120..180` goes
  `+0.56 -> +0.40`.
- Semantic Southern Ocean SON residuals improve from v46 to v48 in the main
  low-pressure target sectors: `-180..-120` goes `+0.94 -> +0.52`,
  `-120..-60` goes `+0.42 -> +0.23`, `0..60` goes `-0.24 -> -0.08`,
  `60..120` goes `+0.19 -> +0.06`, and `120..180` goes `+0.60 -> +0.28`.
- MAM remains controlled after the v48 latitude narrowing.  MAM standardized
  MAE is slightly better than v46, and the wide 45-75S residual does not carry
  the v47 poleward overshoot.
- R2a is still not accepted.  The remaining Southern Ocean blocker is not a
  missing shoulder source anymore; it is residual sector amplitude and
  latitude placement.  In SON, the Pacific/Amundsen and 120-180E sectors remain
  too shallow, while the -60..60 sector still has latitude-dependent sign
  conflict.  Residual North Atlantic compactness remains a secondary watch
  item.

Metrics versus v46:

| Metric | M2 v46 | M2 v48 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.306` | `0.303` |
| standardized-pressure MAE land | `0.416` | `0.415` |
| standardized-pressure MAE ocean | `0.262` | `0.258` |
| standardized-pressure corr all | `0.571` | `0.580` |
| standardized-pressure corr ocean | `0.308` | `0.355` |
| pressure zonal-anomaly corr all | `0.591` | `0.599` |
| pressure zonal-anomaly corr ocean | `0.577` | `0.597` |
| MAM standardized-pressure MAE | `0.355` | `0.354` |
| SON standardized-pressure MAE | `0.430` | `0.417` |
| MAM zonal-anomaly corr | `0.592` | `0.594` |
| SON zonal-anomaly corr | `0.365` | `0.402` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.56s`).

Next M2 step:

1. Keep v48 as the current checkpoint, but do not promote R2a.
2. Continue within Southern Ocean M2 pressure geometry, focusing on SON sector
   amplitude and latitude placement rather than adding broader amplitude.
3. Read Earth SLP, replay pressure, residual, source, transfer, M0 support, and
   M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v49 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v49_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v48:

- Added a SON-only subantarctic north-flank wave transfer term.  This responds
  in the 45-62 S storm-track belt across the southern Atlantic, Pacific, and
  Indian basin margins rather than only inside the semantic Southern Ocean id.
- Added a separate SON-only polar-side trough term in the 58-78 S semantic
  Southern Ocean.  This avoids using one phase for both the north-flank ridge
  and the polar-side trough, which was the main v48 latitude-placement defect.
- Both terms are bounded and tied to M0/M1 support: ocean mask, semantic basin
  exclusion of the Arctic, latitude gates, SST-front support, shelf support,
  open-ocean exposure, and same-latitude SST anomaly where relevant.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- v49 makes the SON Southern Ocean transfer visibly two-layered: a
  subantarctic north-flank ridge/trough wave plus a polar-side trough wave.
- Wide 45-75S SON residuals improve from v48 to v49:
  `-180..-120` goes `+0.465 -> +0.266`, `-60..0` goes
  `-0.174 -> -0.072`, `0..60` goes `-0.158 -> -0.054`,
  `60..120` goes `+0.231 -> +0.160`, and `120..180` goes
  `+0.403 -> +0.123`.
- Semantic Southern Ocean `lat < -55` SON residuals improve from v48 to v49:
  `-180..-120` goes `+0.516 -> +0.349`, `0..60` goes
  `-0.077 -> -0.043`, `60..120` goes `+0.065 -> +0.019`, and
  `120..180` goes `+0.285 -> +0.039`.
- The 45-55S north flank is now much less missing: the west Pacific and
  Australia/New Zealand low sectors deepen, while the South Atlantic /
  South America high sector lifts toward Earth.
- R2a is still not accepted.  Remaining visible owners are:
  high-latitude Southern Ocean trough amplitude in the `-180..0` sectors,
  residual North Atlantic / Nordic / Barents gateway compactness in DJF, and
  subpolar-ocean summer pressure phase that remains too muted.

Metrics versus v48:

| Metric | M2 v48 | M2 v49 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.303` | `0.298` |
| standardized-pressure MAE land | `0.415` | `0.412` |
| standardized-pressure MAE ocean | `0.258` | `0.252` |
| standardized-pressure corr all | `0.580` | `0.593` |
| standardized-pressure corr ocean | `0.355` | `0.402` |
| pressure zonal-anomaly corr all | `0.599` | `0.612` |
| pressure zonal-anomaly corr ocean | `0.597` | `0.625` |
| MAM standardized-pressure MAE | `0.354` | `0.354` |
| SON standardized-pressure MAE | `0.417` | `0.395` |
| MAM zonal-anomaly corr | `0.594` | `0.594` |
| SON zonal-anomaly corr | `0.402` | `0.454` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.13s`).

Next M2 step:

1. Keep v49 as the current checkpoint, but do not promote R2a.
2. Next owner should be selected by map read between:
   high-latitude Southern Ocean `-180..0` trough underdepth,
   DJF North Atlantic / Nordic / Barents gateway compactness, and muted JJA
   subpolar-ocean positive pressure phase.
3. Continue reading Earth SLP, replay pressure, residual, source, transfer,
   M0 support, and M1 support before metrics.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v50 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v50_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v49:

- Added a DJF Atlantic-Arctic gateway transfer term.  It is restricted to the
  Iceland / Greenland / Norwegian / Barents side of the Atlantic-Arctic
  gateway with latitude, longitude, Atlantic/Arctic basin, shelf, and SST-front
  gates.
- Added a small coastal-land inheritance term from the same gateway ocean seed,
  so Iceland/Greenland/Nordic coastal cells can inherit adjacent low-pressure
  support without spreading into broad continental interiors.
- The term explicitly avoids Beaufort and Bering/Chukchi longitudes.  This
  keeps the v45/v49 protection against a spurious pan-Arctic low.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- DJF Greenland Sea residual improves from v49 to v50:
  `+0.566 -> +0.418`.
- DJF Norwegian Sea residual improves:
  `+0.533 -> +0.411`.
- DJF Barents/Kara residual improves:
  `+0.356 -> +0.194`.
- DJF Icelandic core residual improves:
  `+0.283 -> +0.222`; Iceland coastal land improves `+0.467 -> +0.325`.
- DJF Labrador improves modestly `+0.179 -> +0.150`.
- NE Atlantic becomes only slightly more overdeepened `-0.025 -> -0.052`,
  while Beaufort and Bering/Chukchi remain effectively unchanged
  (`-0.095 -> -0.093` and `-0.059 -> -0.057`).
- R2a is still not accepted.  The next visible owner is the missing warm-season
  / shoulder-season subpolar ocean positive-pressure phase.  In MAM and JJA,
  North Pacific / North Atlantic / Arctic subpolar regions show large negative
  residuals while `ocean_pressure_high_source_support` is nearly zero, meaning
  the current M2 high-support logic is too subtropical.

Metrics versus v49:

| Metric | M2 v49 | M2 v50 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.298` | `0.297` |
| standardized-pressure MAE land | `0.412` | `0.412` |
| standardized-pressure MAE ocean | `0.252` | `0.251` |
| standardized-pressure corr all | `0.593` | `0.594` |
| standardized-pressure corr ocean | `0.402` | `0.406` |
| pressure zonal-anomaly corr all | `0.612` | `0.613` |
| pressure zonal-anomaly corr ocean | `0.625` | `0.627` |
| DJF standardized-pressure MAE | `0.229` | `0.227` |
| DJF zonal-anomaly corr | `0.709` | `0.713` |
| SON standardized-pressure MAE | `0.395` | `0.395` |
| SON zonal-anomaly corr | `0.454` | `0.454` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.84s`).

Next M2 step:

1. Keep v50 as the current checkpoint, but do not promote R2a.
2. Next owner should be the warm-season / shoulder-season subpolar ocean
   positive-pressure phase: MAM/JJA North Pacific, North Atlantic, and Arctic
   regions are under-high while the current high-support field is near zero.
3. Preserve the v49 Southern Ocean latitude split and v50 Atlantic-Arctic
   gateway protection against pan-Arctic lows.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v51 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v51_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/replay_ocean_pressure_high_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v50:

- Added MAM/JJA Northern Hemisphere subpolar-ocean high-pressure support.
- MAM support is broad over 45-84 N ocean and shelf/frontal seas, representing
  spring high-latitude ocean/ice-edge pressure support.
- JJA support is narrower and longitude-gated to North Pacific / Gulf of
  Alaska / North Atlantic storm-track margins, intentionally avoiding a
  Beaufort high.
- The support is exported through `ocean_pressure_high_source_support` and
  receives a bounded positive M2 source-to-pressure transfer.  DJF lows,
  Southern Ocean transfer, and Atlantic-Arctic winter gateway terms are left
  unchanged.
- This remains M2-only.  It does not tune R2b wind, R3 currents, SST,
  precipitation, biomes, generated worlds, or downstream acceptance metrics.

Map read:

- The missing MAM/JJA subpolar high-support object class is now visible.
- MAM Aleutian residual improves `-0.291 -> -0.012`; NW Pacific
  `-0.177 -> +0.056`; Gulf Alaska `-0.136 -> +0.110`.
- MAM Icelandic residual improves `-0.865 -> -0.155`; Labrador
  `-0.872 -> -0.228`; NE Atlantic `-0.403 -> +0.084`; Barents
  `-0.590 -> -0.060`; Bering `-0.773 -> -0.047`.
- MAM Greenland Sea and Beaufort remain under-high (`-0.573` and `-0.794`),
  but the missing object class is no longer absent.
- JJA Aleutian residual improves `-0.536 -> -0.452`; Gulf Alaska
  `-0.648 -> -0.558`; Icelandic `-0.439 -> -0.292`; Labrador
  `-0.309 -> -0.159`; NE Atlantic `-0.225 -> -0.092`; Bering
  `-0.214 -> -0.120`.
- JJA Beaufort remains protected (`+0.271 -> +0.271`), so the added high
  support does not create a spurious Arctic summer high.
- R2a is still not accepted.  Remaining visible owners are MAM polar-cap
  latitude/texture regularity, unresolved MAM Greenland/Beaufort under-high,
  JJA North Pacific / Gulf of Alaska under-high, and existing land shoulder-
  season pressure errors.

Metrics versus v50:

| Metric | M2 v50 | M2 v51 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.297` | `0.290` |
| standardized-pressure MAE land | `0.412` | `0.411` |
| standardized-pressure MAE ocean | `0.251` | `0.240` |
| standardized-pressure corr all | `0.594` | `0.614` |
| standardized-pressure corr ocean | `0.406` | `0.499` |
| pressure zonal-anomaly corr all | `0.613` | `0.611` |
| pressure zonal-anomaly corr ocean | `0.627` | `0.621` |
| MAM standardized-pressure MAE | `0.354` | `0.329` |
| MAM zonal-anomaly corr | `0.594` | `0.575` |
| JJA standardized-pressure MAE | `0.212` | `0.207` |
| JJA zonal-anomaly corr | `0.730` | `0.734` |
| DJF standardized-pressure MAE | `0.227` | `0.227` |
| SON standardized-pressure MAE | `0.395` | `0.395` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.71s`).

Next M2 step:

1. Keep v51 as the current checkpoint, but do not promote R2a.
2. Next owner should be selected by map read between MAM polar-cap regularity /
   Greenland-Beaufort under-high, JJA North Pacific / Gulf of Alaska under-high,
   and land shoulder-season pressure errors.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway, and
   v51 Beaufort-protected JJA subpolar high support.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v52 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v52_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_land_pressure_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/replay_m1_energy_support_contact_sheet.png`

Implementation delta from v51:

- Added a Northern Hemisphere high-latitude land shoulder-season phase
  correction in M2 pressure genesis.
- MAM positive support represents residual spring cold / near-freezing shield
  high pressure.  It is gated by high-latitude land, cold seasonal anomaly,
  absolute near-freezing temperature, low interiority, and coast strength.
- SON negative support represents decay of the premature autumn continental
  high caused by the upstream `cooling_tendency` base-pressure term.  It is
  gated by high-latitude land, summer heat memory, autumn cooling memory, low
  elevation, and unfrozen-ground memory.
- The mechanism contributes to both land pressure source and bounded
  source-to-pressure transfer.  It does not modify upstream wind/current/SST
  generation, so R2b/R3/R4 remain frozen as fitting targets.

Map read:

- The root cause is now explicit in code attribution: the upstream base
  pressure proxy contains `-temp_anom + 0.70 * cooling_tendency`.  In MAM,
  warming from winter makes `cooling_tendency` negative and suppresses residual
  continental high pressure.  In SON, cooling from summer makes that term
  strongly positive and creates a premature continental high.  v52 corrects
  only final R2a pressure expression, not the upstream wind-driving base field.
- MAM land residuals improve strongly:
  - Siberia `-0.487 -> -0.114`
  - Eurasia high-latitude land `-0.442 -> +0.025`
  - Canada `-1.145 -> -0.642`
  - North America high-latitude land `-1.013 -> -0.469`
  - Greenland `-0.715 -> +0.047`
- SON premature land highs improve strongly:
  - Siberia `+1.250 -> +0.328`
  - Eurasia high-latitude land `+0.909 -> -0.046`
  - Canada `+1.147 -> +0.109`
  - North America high-latitude land `+1.072 -> +0.177`
  - Greenland `+0.057 -> -0.033`
- Secondary effects are acceptable for this checkpoint: Tibetan MAM residual
  changes slightly (`-0.449 -> -0.422`), Tibetan SON improves
  (`-0.365 -> -0.113`), Sahara MAM/SON both move slightly toward Earth, while
  Australia becomes slightly worse in MAM/SON.
- Visual read: the land support map now has readable MAM/SON land objects, and
  the transfer map carries the intended positive MAM / negative SON high-lat
  land phase.  The MAM polar-cap / Arctic-edge red band remains strong and was
  already present in v51, so v52 is not a promotion point.

Metrics versus v51:

| Metric | M2 v51 | M2 v52 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.290` | `0.284` |
| standardized-pressure MAE land | `0.411` | `0.394` |
| standardized-pressure MAE ocean | `0.240` | `0.239` |
| standardized-pressure corr all | `0.614` | `0.626` |
| standardized-pressure corr land | `0.674` | `0.689` |
| standardized-pressure corr ocean | `0.499` | `0.495` |
| pressure zonal-anomaly corr all | `0.611` | `0.597` |
| pressure zonal-anomaly corr land | `0.604` | `0.590` |
| pressure zonal-anomaly corr ocean | `0.621` | `0.612` |
| MAM standardized-pressure MAE | `0.329` | `0.313` |
| MAM zonal-anomaly corr | `0.575` | `0.566` |
| SON standardized-pressure MAE | `0.395` | `0.388` |
| SON zonal-anomaly corr | `0.454` | `0.436` |
| seasonal wind speed MAE | `1.866 m/s` | `1.866 m/s` |
| wind direction cosine p50 | `0.858` | `0.858` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.73s`).

Next M2 step:

1. Keep v52 as the current checkpoint, but do not promote R2a.
2. Next owner should be MAM polar-cap regularity / Greenland-Beaufort
   under-high or JJA North Pacific / Gulf of Alaska under-high, selected by
   real-Earth map read before code changes.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway,
   v51 Beaufort-protected JJA subpolar high support, and v52 land
   shoulder-phase correction.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v53 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v53_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/replay_ocean_pressure_high_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/replay_pressure_genesis_wave_transfer_seasons.png`

Implementation delta from v52:

- Added a MAM Arctic / Greenland / Beaufort freeze-ocean high-pressure object
  to M2 pressure genesis.
- The support uses near-freezing SST, Arctic/Baffin/Greenland/Beaufort
  longitude gates, shelf support, SST-front support, and high-latitude ocean
  masks.
- Added a Baffin/Labrador gateway sub-support because the Labrador/Baffin
  high-pressure region is lower latitude and was under-triggered by the polar
  cap gate alone.
- The mechanism contributes to ocean high source and bounded source-to-pressure
  transfer.  It leaves wind/current/SST generation unchanged.

Map read:

- MAM Beaufort residual improves `-0.893 -> -0.242`.
- MAM Greenland Sea residual improves `-0.661 -> -0.035`.
- MAM Barents-Kara changes `-0.239 -> +0.105`; this is a mild over-high but
  acceptable because the prior blocker was a much larger Arctic/Greenland/
  Beaufort under-high.
- MAM Baffin/Labrador residual improves `-0.532 -> -0.167`.
- MAM Arctic cap residual improves `-0.662 -> -0.262`.
- Greenland land remains acceptable (`+0.047 -> -0.054`).  Canada MAM land
  remains under-high and should not be hidden by the ocean fix.

Metrics versus v52:

| Metric | M2 v52 | M2 v53 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.284` | `0.281` |
| standardized-pressure MAE land | `0.394` | `0.391` |
| standardized-pressure MAE ocean | `0.239` | `0.236` |
| standardized-pressure corr all | `0.626` | `0.631` |
| standardized-pressure corr ocean | `0.495` | `0.522` |
| pressure zonal-anomaly corr all | `0.597` | `0.601` |
| pressure zonal-anomaly corr ocean | `0.612` | `0.620` |
| MAM standardized-pressure MAE | `0.313` | `0.300` |
| MAM zonal-anomaly corr | `0.566` | `0.587` |
| JJA standardized-pressure MAE | `0.207` | `0.207` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

## M2 Pressure-Genesis v54 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v54_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/replay_ocean_pressure_high_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/replay_pressure_genesis_wave_transfer_seasons.png`

Implementation delta from v53:

- Added a JJA North Pacific / Gulf of Alaska high-pressure object to M2
  pressure genesis.
- The support is restricted to semantic Pacific basin ocean, 41-67 N, and Gulf
  Alaska / Aleutian longitude gates, with shelf, SST-front, and cool
  same-latitude SST support.
- The support is zero over Beaufort/Arctic so v51's Beaufort-protected summer
  behavior remains intact.

Map read:

- JJA Gulf Alaska residual improves `-0.585 -> -0.268`.
- JJA Aleutian residual improves `-0.446 -> -0.156`.
- JJA North Pacific residual improves `-0.376 -> -0.121`.
- JJA NW Pacific changes `-0.117 -> +0.139`, a mild over-high tradeoff.
- JJA Beaufort remains protected (`+0.308 -> +0.303`).
- MAM and SON metrics remain unchanged from v53, as intended.
- R2a is still not accepted.  The next visible owner should be selected by
  map read between MAM Canada land under-high, residual MAM Arctic cap
  under-high, and SON / Southern Ocean pressure-wave residuals.

Metrics versus v53:

| Metric | M2 v53 | M2 v54 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.281` | `0.279` |
| standardized-pressure MAE land | `0.391` | `0.391` |
| standardized-pressure MAE ocean | `0.236` | `0.234` |
| standardized-pressure corr all | `0.631` | `0.633` |
| standardized-pressure corr ocean | `0.522` | `0.530` |
| pressure zonal-anomaly corr all | `0.601` | `0.604` |
| pressure zonal-anomaly corr ocean | `0.620` | `0.624` |
| MAM standardized-pressure MAE | `0.300` | `0.300` |
| JJA standardized-pressure MAE | `0.207` | `0.203` |
| JJA zonal-anomaly corr | `0.734` | `0.741` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.66s`).

Next M2 step:

1. Keep v54 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: MAM Canada land under-high,
   residual MAM Arctic cap under-high, or SON / Southern Ocean pressure-wave
   residuals.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway,
   v51 Beaufort-protected JJA subpolar high support, v52 land shoulder-phase
   correction, v53 MAM Arctic freeze-ocean high, and v54 JJA North Pacific
   high support.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v55 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v55_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/replay_land_pressure_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/replay_pressure_genesis_wave_transfer_seasons.png`

Implementation delta from v54:

- Added a MAM North America spring land-high object in M2 pressure genesis.
- The support is restricted to 43-73 N and -165..-45 longitude, with a compact
  Canada-centered longitude gate, low-elevation support, cold / near-freezing
  memory, low-interiority weighting, and coast-strength modulation.
- Rejected generic candidate: a low-interiority / low-elevation high-latitude
  land support improved Canada but also raised Eurasia/Siberia, worsening MAM
  land MAE.  v55 therefore uses a real-Earth regional gate for this replay
  checkpoint rather than global coefficient tuning.
- Fixed the land-support diagnostic merge to use `np.maximum.reduce`, removing
  the NumPy warning caused by a three-argument `np.maximum` call.

Map read:

- MAM Canada residual improves `-0.665 -> -0.321`.
- MAM North America high-latitude land residual improves `-0.502 -> -0.243`.
- MAM Greenland land remains near target (`-0.054 -> -0.009`).
- MAM Siberia and Eurasia remain near target (`-0.107 -> -0.123`,
  `+0.008 -> -0.009`), so the fix does not reintroduce a broad Eurasian
  over-high.
- MAM Alaska becomes mildly over-high (`+0.066 -> +0.244`), retained as an
  acceptable tradeoff because Canada/North America and MAM metrics improve.
- MAM Arctic ocean changes slightly (`-0.262 -> -0.281`), so residual Arctic
  cap under-high remains a possible later owner.
- JJA and SON pressure metrics remain unchanged from v54, as intended.

Metrics versus v54:

| Metric | M2 v54 | M2 v55 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.279` | `0.278` |
| standardized-pressure MAE land | `0.391` | `0.387` |
| standardized-pressure MAE ocean | `0.234` | `0.234` |
| standardized-pressure corr all | `0.633` | `0.635` |
| standardized-pressure corr land | `0.684` | `0.688` |
| standardized-pressure corr ocean | `0.530` | `0.529` |
| pressure zonal-anomaly corr all | `0.604` | `0.607` |
| pressure zonal-anomaly corr land | `0.593` | `0.600` |
| pressure zonal-anomaly corr ocean | `0.624` | `0.623` |
| MAM standardized-pressure MAE | `0.300` | `0.296` |
| MAM zonal-anomaly corr | `0.587` | `0.604` |
| JJA standardized-pressure MAE | `0.203` | `0.203` |
| SON standardized-pressure MAE | `0.388` | `0.388` |

Tests:

- Targeted regression tests pass without the previous NumPy warning:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.57s`).

Next M2 step:

1. Keep v55 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: residual MAM Arctic cap
   under-high or SON / Southern Ocean pressure-wave residuals.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway,
   v51 Beaufort-protected JJA subpolar high support, v52 land shoulder-phase
   correction, v53 MAM Arctic freeze-ocean high, v54 JJA North Pacific high,
   and v55 North America spring land-high support.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v62 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v62_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/pressure_standardized_delta_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/replay_land_pressure_source_support_seasons.png`

Implementation delta from v60:

- Added SON boreal autumn land-high relief support.  It is controlled by
  summer heat memory, autumn cooling memory, low elevation, unfrozen-ground
  memory, coast / low-interiority leakage, and terrain-barrier escape.  The
  intent is to prevent broad cold-high pressure from appearing too early over
  low-elevation high-latitude land.
- Added SON North Atlantic / Icelandic autumn low support.  It is restricted
  to semantic Atlantic ocean cells, uses shelf / SST-front / same-latitude SST
  support, and tapers westward to avoid turning Labrador into a false low.
- Added a SON Southern Ocean sector-low wave adjustment after v61 showed that
  the land / North Atlantic fix raised Southern Ocean residuals.  The v62
  adjustment is strongest in Pac/Amundsen and Indian sectors and cuts out most
  Atlantic-sector support.

Map read:

- SON Siberia residual improves `+0.338 -> +0.137`.
- SON Canada improves `+0.204 -> -0.057`.
- SON Alaska improves `+0.395 -> +0.078`.
- SON Greenland moves `-0.123 -> -0.196`; this is a mild cost but remains
  smaller than the removed Siberia / Canada / Alaska over-high.
- SON North Atlantic subpolar improves `+0.247 -> +0.074`.
- SON Icelandic sector improves `+0.455 -> +0.098`.
- SON Labrador is protected (`+0.076 -> +0.078`).
- SON North Pacific subpolar improves `-0.117 -> -0.045`.
- Southern Ocean all is preserved at the v60 residual level
  (`+0.066 -> +0.066`) after the v61 side effect.
- Southern Ocean Pac/Amundsen is near v60 (`+0.106 -> +0.109`), Atlantic is
  near target (`-0.026 -> -0.013`), Indian improves (`+0.183 -> +0.133`), and
  Aus-Pac changes mildly (`-0.015 -> +0.020`).
- MAM central Arctic / Baffin-Labrador and JJA North Atlantic / NW Pacific
  were intentionally unchanged and remain later owners.

Metrics versus v60:

| Metric | M2 v60 | M2 v62 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.2748` | `0.2723` |
| standardized-pressure MAE land | `0.3809` | `0.3770` |
| standardized-pressure MAE ocean | `0.2320` | `0.2301` |
| standardized-pressure corr all | `0.6417` | `0.6441` |
| standardized-pressure corr land | `0.691` | `0.693` |
| standardized-pressure corr ocean | `0.543` | `0.547` |
| pressure zonal-anomaly corr all | `0.615` | `0.611` |
| pressure zonal-anomaly corr ocean | `0.631` | `0.624` |
| SON standardized-pressure MAE | `0.3878` | `0.3779` |
| SON zonal-anomaly corr | `0.436` | `0.416` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.30s`).

Next M2 step:

1. Keep v62 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: MAM Baffin-Labrador /
   central Arctic under-high, JJA North Atlantic subpolar under-high / NW
   Pacific over-high, or residual SON high-latitude texture.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway,
   v51 Beaufort-protected JJA subpolar high support, v52 land shoulder-phase
   correction, v53 MAM Arctic freeze-ocean high, v54/v60 North Pacific high
   support, v55/v59 North America / Canadian spring highs, v58 DJF lows /
   North America winter relief, and v62 SON land / North Atlantic / Southern
   Ocean corrections.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.

## M2 Pressure-Genesis v63 Checkpoint

Date: 2026-07-07

Evidence packet:

- `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v63_20260707/`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/real_earth_pressure_replay_contact_sheet.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/pressure_standardized_delta_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/replay_pressure_genesis_source_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/replay_pressure_genesis_wave_transfer_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/replay_ocean_pressure_low_source_support_seasons.png`
- `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/replay_land_pressure_source_support_seasons.png`

Implementation delta from v62:

- Added a JJA Eurasian summer thermal-low object.  The support is restricted
  to the heated low-elevation Eurasian belt from Arabia / Iran through India
  and East Asia.  It is gated by JJA thermal anomaly, low elevation,
  continent / coast expression, and terrain-barrier escape.  North America is
  excluded because v62 already over-deepens its JJA land low.
- Added a JJA western Pacific marginal-sea low / trough object.  The support
  is restricted to semantic Pacific / Indian marginal-sea cells in the
  118E-178E, 24-64N domain and gated by shelf, SST-front, and same-latitude
  SST support.

Map read:

- JJA NW Pacific residual improves `+0.360 -> +0.070`.
- JJA Kuroshio / Oyashio improves `+0.312 -> +0.065`.
- JJA Japan / East China Sea improves `+0.747 -> +0.606`; this remains
  visibly too high and is a possible later owner.
- JJA East Asia land improves `+0.386 -> +0.225`.
- JJA NE Asia land improves `+0.223 -> +0.072`.
- JJA China lowland improves `+0.539 -> +0.352`.
- JJA India improves `+0.431 -> +0.312`.
- JJA Arabia / Iran improves `+0.636 -> +0.499`, still underexpressed.
- JJA North America west improves mildly as a normalization side effect
  (`-0.568 -> -0.515`) despite receiving no direct support.
- JJA Gulf Alaska is preserved (`-0.132 -> -0.128`), central North Pacific is
  preserved (`-0.172 -> -0.153`), and MAM central Arctic / SON target regions
  are unchanged.

Metrics versus v62:

| Metric | M2 v62 | M2 v63 |
| --- | ---: | ---: |
| standardized-pressure MAE all | `0.2723` | `0.2705` |
| standardized-pressure MAE land | `0.3770` | `0.3728` |
| standardized-pressure MAE ocean | `0.2301` | `0.2293` |
| standardized-pressure corr all | `0.6441` | `0.6493` |
| standardized-pressure corr land | `0.6929` | `0.6989` |
| standardized-pressure corr ocean | `0.5469` | `0.5510` |
| pressure zonal-anomaly corr all | `0.611` | `0.618` |
| pressure zonal-anomaly corr ocean | `0.624` | `0.630` |
| JJA standardized-pressure MAE | `0.1995` | `0.1923` |
| JJA zonal-anomaly corr | `0.749` | `0.778` |

Tests:

- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 105.55s`).

Next M2 step:

1. Keep v63 as the current checkpoint, but do not promote R2a.
2. Select the next owner by real-Earth map read: MAM Baffin-Labrador /
   central Arctic under-high, JJA North Atlantic / Icelandic under-high, JJA
   Japan / East China Sea and Arabia / Iran summer low underexpression, or
   residual SON high-latitude texture.
3. Preserve v49 Southern Ocean latitude split, v50 Atlantic-Arctic gateway,
   v51 Beaufort-protected JJA subpolar high support, v52 land shoulder-phase
   correction, v53 MAM Arctic freeze-ocean high, v54/v60 North Pacific high
   support, v55/v59 North America / Canadian spring highs, v58 DJF lows /
   North America winter relief, v62 SON land / North Atlantic / Southern
   Ocean corrections, and v63 JJA Eurasian / western Pacific corrections.
4. R2b wind, R3 currents, R4 SST, precipitation, biomes, and generated worlds
   remain blocked until R2a pressure passes map-read acceptance.
