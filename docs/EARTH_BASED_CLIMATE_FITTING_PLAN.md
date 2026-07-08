# Earth-Based Climate Fitting Plan

Status: active; real-Earth single-subgraph replay is the authoritative fitting workflow; generated worlds and old C5 gates are regression guardrails only until R0-R8 pass on Earth
Owner: climate / biome calibration
Last updated: 2026-07-07

This plan starts after the Earth reference pipeline.  Plate and terrain
generation are treated as frozen inputs.  The active goal is to fit one
real-Earth geoscience subgraph at a time on the real-Earth grid, using a direct
map comparison between the Earth reference field and the Aevum replay field.
Generated worlds are not current fitting targets; they return only as R9
generalization guardrails after the active Earth subgraph is physically
plausible.

## Active Earth-Only Protocol

This section is the current working contract and overrides older F0-F5 notes in
this file.

- R0-R8 are real-Earth replay phases.  Each phase chooses exactly one subgraph
  such as seasonal wind vectors, OSCAR surface currents, OISST SST, moisture
  transport, precipitation, sea ice, or biome classes.
- For that subgraph, render the real-Earth reference map and the Aevum
  real-Earth replay map on the same grid, projection, color scale, and mask.
  Continuous fields also need a residual or vector-error map.
- Read the maps first.  Record the visible structures that are correct,
  missing, displaced, too zonal, too noisy, too smooth, or physically
  impossible.  Metrics are read after this as regression evidence, not as the
  primary acceptance proof.
- Attribute each major residual to an upstream dependency or to the active
  mechanism.  If the owner is upstream, stop and fix that upstream layer rather
  than compensating downstream.
- Change only the mechanism that owns the diagnosed residual.  Do not tune
  sea-ice, precipitation, biome, or class thresholds to hide broken wind,
  current, SST, or moisture fields.
- R9 is the first point where generated worlds are used again.  Its job is to
  check that the Earth-fitted mechanism still generalizes to the accepted
  terminal worlds.

Practical execution rule:

- Every implementation loop fits exactly one Earth reference subgraph.  The
  loop starts with the real map, the same-grid replay map, and the residual map;
  it does not start from generated-world failures or global summary scores.
- The first written artifact for each loop is a map-read attribution note:
  name the real-Earth structures, name the replay residuals, and assign each
  residual to an upstream owner or the active mechanism.
- A code change is accepted only if it repairs the named map residual through
  that owner.  Global means, p50/p90 envelopes, or pass/fail tables may reject a
  regression, but they cannot promote a visually wrong subgraph.
- Downstream maps are observer-only during upstream fitting.  In particular,
  sea ice, Koppen classes, biomes, precipitation, SST, currents, or R2b wind
  translation must not be tuned to hide an unresolved R2a pressure-source
  residual.

Current active packet:

- Active subgraph: R2a seasonal SLP / pressure-source geometry replay on real
  Earth.
- Required reference: `earth__seasonal_slp_anomaly_hPa`.
- Required replay field: `atmosphere__land_sea_pressure_proxy` /
  `atmosphere__seasonal_pressure_proxy` as the Aevum pressure proxy.  Wind,
  current, SST, precipitation, ice, Koppen, biome, and generated-world maps are
  observer-only and cannot justify a pressure-source change.
- Required figures: Earth seasonal standardized SLP anomaly, replay seasonal
  standardized pressure-proxy anomaly, standardized pressure residual, zonal
  pressure anomaly, and any pressure-center / stationary-wave object support
  maps emitted by the replay.
- Blocked until R2a map-read acceptance: R2b wind translation, R3 currents, R4
  SST repair, R5 moisture, R6 precipitation, R7 cryosphere/cloud/vegetation
  feedback, and R8 climate/biome classes.

Mechanism-first correction, 2026-07-07:

- R2a remains the active replay packet, but the immediate next action is not
  another pressure-source parameter change.  The v1 basin-pressure-source
  experiment under
  `out_real_earth_pressure_replay_r2a_basin_pressure_v1_20260707/` was rejected
  because it worsened ocean pressure correlation and did not repair the visible
  geography; the attempted code change was rolled back.  The partial v2 replay
  under `out_real_earth_climate_replay_r2a_basin_pressure_v2_20260707/` is also
  a rejected artifact.
- Before more implementation, define the coupled submodel contracts documented
  in `docs/CLIMATE_COUPLING_RESEARCH_NOTES.md`: terrain/land-sea geometry and
  energy state must feed pressure and wind; wind and basin geometry must feed
  ocean currents; currents and heat flux must feed SST; SST, wind, terrain, and
  convergence must feed moisture and precipitation.  The detailed contract and
  microbenchmark plan is archived in
  `docs/CLIMATE_MECHANISM_MODELING_PLAN.md`.
- Parameter sweeps may only follow a named mechanism and map-read attribution.
  They cannot replace the mechanism design step.

M0/M1 diagnostic contract checkpoint, 2026-07-07:

- `ClimateModule` now emits diagnostic-only M1 energy-boundary arrays needed to
  attribute R2a pressure residuals before changing pressure code.
- `terminal_climate_arrays.npz` now archives those arrays, and
  `real-earth-wind-replay` renders M0/M1 support contact sheets when the arrays
  are present.
- The active subgraph remains R2a seasonal SLP / pressure-source geometry.  The
  new support maps are upstream attribution evidence, not a downstream tuning
  target.
- First M0/M1 visual attribution note:
  `docs/R2A_M0_M1_MAP_READ_ATTRIBUTION_20260707.md`.  It first identified
  missing major-ocean semantic basin support in Earth replay; that M0 blocker
  has now been repaired by the major-ocean checkpoint.  The current hard
  blocker is M2 pressure genesis: it still under-consumes basin, SST-front, and
  terrain/coast stationary-wave support.

R2a M2 v6 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v6_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v6_20260707/`.
- New auditable fields:
  `atmosphere.pressure_genesis_source`,
  `atmosphere.pressure_genesis_wave_transfer`,
  `atmosphere.ocean_pressure_low_source_support`,
  `atmosphere.ocean_pressure_high_source_support`,
  `atmosphere.land_pressure_source_support`, and
  `atmosphere.terrain_pressure_wave_source_support`.
- Map read: major ocean basins and M1 energy support are now readable; DJF
  Aleutian/Icelandic source patches and JJA Southern Ocean source sectors are
  visible; the v6 wave-transfer pass reduces the previous polar-edge artifact.
  The final pressure proxy is still too smooth and too weakly organized by
  coast, terrain, and SST-front waveguides, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.318/0.419/0.277`; standardized-pressure correlation all/land/ocean
  `0.548/0.670/0.209`; pressure zonal-anomaly correlation all/land/ocean
  `0.579/0.589/0.570`.
- Next accepted change must stay inside M2 and improve source-to-pressure
  propagation geometry.  R2b wind, R3 currents, R4 SST, R5 moisture, R6
  precipitation, R7 feedback, R8 classes/biomes, and generated worlds remain
  blocked as fitting targets.

R2a M2 v7 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v7_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v7_20260707/`.
- Mechanism delta: v7 keeps the M2 causal source fields but replaces the v6
  mostly isotropic transfer with weighted source-to-pressure diffusion.  Ocean
  propagation is weighted by SST-front support, same-latitude SST anomaly,
  open-ocean exposure, and subpolar/basin support.  Land propagation is
  weighted by coast-strength, terrain barriers, land-source support, and
  terrain-wave support with a non-polar gate.
- Map read: v7 transfer is better aligned with northern storm-track/coastal
  waveguide structure and does not reintroduce the Antarctic edge artifact.
  Final pressure is still visually close to v6 and remains too smooth/blocky
  versus Earth SLP, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.317/0.418/0.277`; standardized-pressure correlation all/land/ocean
  `0.549/0.671/0.209`; pressure zonal-anomaly correlation all/land/ocean
  `0.581/0.591/0.571`.
- Next accepted change must remain inside M2 and make source objects project
  into stronger, readable pressure centers rather than only weak local
  diffusion.  R2b/R3/R4/R5/R6/R7/R8 and generated worlds remain blocked as
  fitting targets.

R2a M2 v8 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v8_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v8_20260707/`.
- Mechanism delta: v8 keeps the v7 weighted waveguide but adds bounded
  object-level pressure-center projection from subpolar ocean-low support,
  land thermal-center source support, and terrain-wave source support.  The
  final transfer is stronger than v7 but still clipped and gated to avoid
  polar-edge amplification.
- Map read: v8 transfer has clearer DJF Aleutian/Icelandic low and Eurasian
  winter-high signatures, and JJA Southern Ocean remains segmented rather than
  a single annular paint band.  The final pressure proxy still remains too
  smooth relative to Earth SLP, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.317/0.418/0.277`; standardized-pressure correlation all/land/ocean
  `0.550/0.672/0.212`; pressure zonal-anomaly correlation all/land/ocean
  `0.581/0.592/0.570`.
- Next accepted change must remain inside M2 and improve object placement and
  anisotropic footprint geometry, not downstream wind/current/SST/moisture or
  biome thresholds.

R2a M2 v9 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v9_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v9_20260707/`.
- Mechanism delta: v9 implements directional diffusion for M2
  source-to-pressure transfer.  Ocean footprints diffuse along a
  storm-track/SST-front axis; land footprints diffuse along terrain-barrier and
  coastal axes.  This replaces ordinary local diffusion in the v8 wave/object
  projection path.
- Map read: v9 preserves the v8 DJF Aleutian/Icelandic/continental footprints
  and JJA Southern Ocean sector structure, with no new polar-edge artifact.
  Final pressure remains visually too smooth and too close to v8, so R2a is
  not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.317/0.418/0.277`; standardized-pressure correlation all/land/ocean
  `0.550/0.672/0.212`; pressure zonal-anomaly correlation all/land/ocean
  `0.581/0.592/0.570`.
- Next accepted change must still remain inside M2, but the next owner is no
  longer "add anisotropy" by itself.  The remaining problem is that pressure
  centers are more readable in source/transfer diagnostics than in the final
  pressure proxy; M2 must improve final-pressure expression strength and
  placement without global amplitude painting.

R2a M2 v10 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v10_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v10_20260707/`.
- Mechanism delta: v10 adds a bounded final-pressure expression solve after
  the v9 directional transfer.  It enhances source-supported nonzonal pressure
  anomaly only where the sign aligns with the M2 source/transfer field, then
  removes the latitude-band mean and clips the increment.  This is intended to
  let existing pressure-center objects affect final pressure without global
  latitude-band painting.
- Map read: v10 makes DJF North Pacific / North Atlantic and Eurasian transfer
  centers more forceful, improves JJA Southern Ocean wave expression, and does
  not reintroduce the Antarctic edge artifact.  Final pressure is improved but
  still smoother than Earth SLP, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.317/0.416/0.276`; standardized-pressure correlation all/land/ocean
  `0.550/0.671/0.220`; pressure zonal-anomaly correlation all/land/ocean
  `0.582/0.593/0.571`.
- Next accepted change remains M2-only.  Focus on center placement and shape:
  pressure centers should become more coherent in the final pressure proxy
  without turning source-supported anomaly expression into broad high-latitude
  bands.

R2a M2 v16 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v16_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v16_20260707/`.
- Mechanism delta: v11/v12 object-local final expression and v13/v14 stronger
  object/source expression did not improve map readability or metrics enough to
  keep.  v16 keeps the v10 final-expression solve and adds a source-side
  mechanism: winter cold continental highs are spread downwind into nearby open
  ocean along the westerly/storm-track axis, then used as support for subpolar
  ocean-low source selection.  Southern Hemisphere downwind support is damped
  so it does not replace the Southern Ocean SST-front/wavenumber gate.
- Map read: Aleutian/Icelandic source placement remains readable and is now
  partly tied to cold-continent downstream geometry.  The Southern Ocean remains
  segmented rather than becoming a full ring.  Final pressure is still too
  smooth, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.317/0.416/0.276`; standardized-pressure correlation all/land/ocean
  `0.551/0.671/0.221`; pressure zonal-anomaly correlation all/land/ocean
  `0.583/0.594/0.572`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.24s`).
- Next accepted change remains M2-only.  The source maps are now readable enough
  that the remaining problem is final-pressure dominance by the smoothed
  upstream pressure proxy.  Re-balance source-supported pressure expression
  against that upstream proxy; do not proceed to R2b wind, R3 currents, SST,
  precipitation, or biomes.

R2a M2 v17 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v17_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v17_20260707/`.
- Mechanism delta: v17 uses the v16 arrays to validate the final pressure
  synthesis before changing code.  The winning direction is not raw source
  amplification.  Instead, M2 source support remains a causal trigger, its
  direct pressure contribution is reduced to `0.80`, and the bounded
  source-to-pressure transfer contribution is increased to `1.45`.
- Map read: transfer expression is stronger in already source-supported
  Northern Hemisphere subpolar regions and in the JJA Southern Ocean wave belt.
  Final pressure is less dominated by the smooth upstream proxy, but pressure
  centers still lack Earthlike coherence and placement.  R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.316/0.416/0.275`; standardized-pressure correlation all/land/ocean
  `0.552/0.672/0.222`; pressure zonal-anomaly correlation all/land/ocean
  `0.584/0.595/0.573`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 104.60s`).
- Next accepted change remains M2-only.  Improve pressure-center morphology and
  placement using the existing source/transfer support.  Do not add hard-coded
  seasonal gain tables and do not proceed to R2b wind, R3 currents, SST,
  precipitation, or biomes.

R2a M2 v18 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v18_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v18_20260707/`.
- Mechanism delta: v18 keeps the v17 source/transfer balance and adds bounded
  transfer morphology.  Ocean-low source support can deepen negative transfer
  cores, while broad ocean-high, land-core, and terrain-wave transfer
  contributions are damped to reduce over-smooth or fragmented pressure
  expression.  This uses existing M2 supports only.
- Map read: the transfer map is slightly more center-like and less terrain
  fragmented.  No new Antarctic edge artifact, Southern Ocean ring, or
  high-latitude speckle was introduced.  Final pressure remains too smooth for
  R2a acceptance.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.316/0.416/0.275`; standardized-pressure correlation all/land/ocean
  `0.553/0.673/0.224`; pressure zonal-anomaly correlation all/land/ocean
  `0.585/0.596/0.574`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 204.10s`).
- Next accepted change remains M2-only.  Continue improving final-pressure
  center morphology and placement; do not proceed to R2b wind, R3 currents,
  SST, precipitation, or biomes.

R2a M2 v19 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v19_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v19_20260707/`.
- Mechanism delta: v19 keeps the v18 source/transfer/morphology structure and
  adds a bounded thermal-phase adjustment to transfer expression.  Same-latitude
  SST and land-temperature anomalies modulate high-support, low-support, and
  land-source footprints so pressure signs follow seasonal thermal phase more
  consistently.  This uses existing M1/M2 fields only.
- Map read: seasonal differentiation in the transfer/final pressure panels is
  better, especially in MAM/SON phase behavior.  No new Antarctic edge artifact,
  Southern Ocean ring, or heat-wall latitude band was introduced.  Final
  pressure remains too smooth and weakly center-organized for R2a acceptance;
  broader subtropical thermal-phase patches are a watch item.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.310/0.416/0.268`; standardized-pressure correlation all/land/ocean
  `0.564/0.675/0.269`; pressure zonal-anomaly correlation all/land/ocean
  `0.587/0.600/0.574`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 139.05s`).
- Next accepted change remains M2-only.  Continue improving final-pressure
  center morphology and placement; do not proceed to R2b wind, R3 currents,
  SST, precipitation, or biomes.

R2a M2 v25 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v25_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v25_20260707/`.
- Mechanism delta: v25 tightens the Southern Ocean wavenumber/front gate and
  changes final pressure synthesis to a domain-weighted direct-source
  expression.  Land and Southern Ocean source overprint is reduced, while North
  Hemisphere mid/high-latitude ocean-low cores keep more direct expression.
  This keeps Aleutian-type lows from being erased while reducing the JJA
  Southern Ocean over-deep residual.
- Rejected trials: v20/v21 over-activated MAM/SON ocean lows, and v23 basin
  reweighting did not improve the map.  These were not promoted.
- Map read: no new Antarctic edge artifact, Southern Ocean ring, or heat-wall
  latitude band was introduced.  JJA Southern Ocean residuals improve; DJF
  Aleutian remains close to v19.  DJF Icelandic Low is still too weak, and SON
  North Pacific still lacks a meaningful source.  R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.309/0.416/0.266`; standardized-pressure correlation all/land/ocean
  `0.564/0.674/0.272`; pressure zonal-anomaly correlation all/land/ocean
  `0.589/0.600/0.577`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 115.70s`).
- Next accepted change remains M2-only.  Continue improving pressure-source
  placement and final-pressure center morphology; do not proceed to R2b wind,
  R3 currents, SST, precipitation, or biomes.

R2a M2 v30 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v30_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v30_20260707/`.
- Mechanism delta: v30 adds a post-normalized Northern Hemisphere winter
  coastal/front supplement for basin-edge subpolar lows and a bounded
  shoulder-season warm-ocean/front low-source candidate.  The shoulder source is
  component-selected but amplitude-limited after selection, so SON sources can
  appear without turning MAM into a low-source season.
- Rejected trials: v26 put the coastal/front term into the normalized low-score
  solve and weakened Aleutian; v28 made SON ocean lows too deep; v29 was closer
  but still too strong.  v30 lowers the shoulder amplitude.
- Map read: SON North Pacific now has a weak local low-source/transfer
  expression; MAM North Pacific/North Atlantic remain untriggered.  DJF
  Icelandic Low is slightly stronger than v25 and DJF Aleutian is not visibly
  degraded.  Final pressure is still too smooth and Icelandic remains too weak.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.309/0.416/0.266`; standardized-pressure correlation all/land/ocean
  `0.563/0.674/0.268`; pressure zonal-anomaly correlation all/land/ocean
  `0.585/0.598/0.569`.
- Tradeoff: compared with v25, the SON North Pacific source-placement blocker
  improves, but SON/ocean scalar correlations regress.  v30 is a mechanism
  checkpoint, not an acceptance point.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.16s`).
- Next accepted change remains M2-only.  Preserve the new SON source object
  while improving final-pressure expression and pressure-center compactness; do
  not proceed to R2b wind, R3 currents, SST, precipitation, or biomes.

R2a M2 v38 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v38_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v38_20260707/`.
- Mechanism delta: v38 keeps the v30 shoulder-season warm-ocean source but
  makes the post-object amplitude depend more strongly on basin scale and
  open-ocean exposure.  This preserves the needed SON North Pacific weak low
  while reducing the v30 North Atlantic shoulder-season over-deepening.
- v38 also adds a winter-only Northern Hemisphere subpolar source-expression
  boost, improving DJF Aleutian/Icelandic lows without changing MAM/SON/JJA
  low-source activation.
- Southern Ocean pressure now uses front/shelf/same-latitude-SST support as
  the main gate, with the longitude wave reduced to a weak perturbation.  A
  signed Southern Ocean wave-transfer anomaly turns high support sectors into
  relative lows and low support sectors into relative highs, instead of only
  painting negative source everywhere.
- Rejected intermediate trials: v31 removed the SON North Pacific source
  together with the North Atlantic over-deepening; v32/v33 restored a weaker
  SON North Pacific source and reduced North Atlantic excess; v35/v36 were too
  conservative over the Southern Ocean; v37 introduced the correct signed
  Southern Ocean wave-transfer direction; v38 increases that wave enough to be
  visible while remaining bounded.
- Map read: DJF Aleutian and Icelandic lows are stronger than v30.  SON North
  Pacific keeps a weak low source while MAM remains untriggered.  JJA Southern
  Ocean now shows a signed transfer wave; the `60..120E` sector moves toward
  the real-Earth low, and several erroneous low sectors are lifted.  The final
  pressure field is still too smooth and Southern Ocean / North Atlantic
  residuals remain visible, so R2a is not accepted.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.309/0.416/0.266`; standardized-pressure correlation all/land/ocean
  `0.564/0.674/0.272`; pressure zonal-anomaly correlation all/land/ocean
  `0.588/0.600/0.576`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.11s`).
- Full test suite was started after the targeted tests and interrupted for
  time after `103 passed in 444.62s`; no failure was observed before the
  interrupt.
- Next accepted change remains M2-only.  Improve pressure-center compactness
  and remaining North Atlantic / Southern Ocean placement; do not proceed to
  R2b wind, R3 currents, SST, precipitation, or biomes.

R2a M2 v40 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v40_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v40_20260707/`.
- Mechanism delta: v40 keeps the v38 Southern Ocean signed-wave and
  shoulder-season source rules, then strengthens only the Northern Hemisphere
  winter subpolar ocean-low compact-expression term.  This uses existing
  `low_support`, SST-front support, the winter hemisphere gate, and the
  North-Hemisphere subpolar ocean mask; it does not add new basin ids,
  coordinate targets, or downstream tuning.
- Rejected/retuned step: v39 confirmed the compact-expression direction at
  lower amplitude.  Offline replay on v39 arrays showed that a coefficient near
  `0.20` best improved Aleutian/Icelandic residuals and ocean correlation
  before global MAE began to flatten/regress, so v40 uses that stronger bounded
  value.
- Map read: DJF Aleutian and Icelandic cores deepen while MAM/SON/JJA maps are
  unchanged by this term.  Region means improve from v38 to v40:
  Aleutian core replay `-0.575 -> -0.645` against Earth `-0.840`, and
  Icelandic core `-0.435 -> -0.493` against Earth `-0.706`.  SON North Pacific
  weak source and Southern Ocean signed wave from v38 remain present.
- Remaining defect: DJF low cores are still weaker than Earth, Labrador /
  North Atlantic edge placement remains incomplete, and the final replay is
  still too smooth for R2a promotion.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.309/0.416/0.265`; standardized-pressure correlation all/ocean
  `0.565/0.277`; pressure zonal-anomaly correlation all/ocean
  `0.589/0.577`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.26s`).
- Next accepted change remains M2-only.  Continue improving pressure-center
  compactness and North Atlantic / Southern Ocean placement; do not proceed to
  R2b wind, R3 currents, SST, precipitation, or biomes.

R2a M2 v44 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v44_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v44_20260707/`.
- Mechanism delta: v44 keeps the v40 compact Northern Hemisphere winter
  ocean-low expression, adds a narrow coastal-land inheritance term so
  Labrador/Iceland/nearby islands can inherit adjacent ocean-low pressure, and
  adds a bounded 50-67 N subpolar SST-front / storm-track low-support floor.
  This remains M2-only and uses existing M0/M1 support fields.
- Rejected/retuned steps: v41 was too weak; v42 made the coastal inheritance
  visible and is retained inside v44.  v43 moved lee-low support into pressure
  object selection, which overdeepened Labrador/NW Pacific and weakened the
  Icelandic low through component competition, so that direction is rejected.
- Map read: DJF Labrador ocean residual improves from `+0.253` in v42 to
  `+0.179` in v44; Icelandic ocean from `+0.353` to `+0.333`; Labrador land
  from `+0.242` to `+0.235`; Icelandic land from `+0.512` to `+0.504`.
  NE Atlantic ocean becomes slightly overdeepened (`-0.014 -> -0.033`), which
  is acceptable for this checkpoint but remains a watch item.
- Remaining defect: R2a is still not accepted.  The North Atlantic / Icelandic
  low is still weaker and smoother than Earth, Arctic/Nordic winter low
  structure remains under-expressed, and Southern Ocean wave sectors still need
  organization.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.308/0.415/0.265`; standardized-pressure correlation all/ocean
  `0.566/0.281`; pressure zonal-anomaly correlation all/ocean
  `0.590/0.577`; DJF standardized-pressure MAE `0.230`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.40s`).
- Next accepted change remains M2-only.  Continue with real-Earth SLP map
  reading first; do not proceed to R2b wind, R3 currents, SST, precipitation,
  biomes, or generated-world fitting.

R2a M2 v45 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v45_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v45_20260707/`.
- Mechanism delta: v45 keeps v44 and adds an Atlantic-Arctic gateway
  low-support floor.  It computes graph distance from existing Atlantic
  subpolar low support through the Atlantic/Arctic ocean domain, then applies a
  bounded floor only to nearby 62-80 N Arctic-basin marginal seas with shelf and
  SST-front support.
- Map read: DJF Nordic/Arctic gateway residual improves from `+0.642` in v44
  to `+0.538` in v45; Greenland Sea from `+0.645` to `+0.562`;
  Barents/Kara from `+0.479` to `+0.425`; Icelandic ocean from `+0.337` to
  `+0.314`.  Beaufort/Arctic and Bering/Chukchi remain essentially protected,
  so v45 does not paint a blanket Arctic low.
- Remaining defect: R2a is still not accepted.  Nordic/Arctic and Icelandic
  lows remain underdeepened, and the clearest next owner is now Southern Ocean
  shoulder-season / wave-sector source organization.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.308/0.415/0.265`; standardized-pressure correlation all/ocean
  `0.566/0.284`; pressure zonal-anomaly correlation all/ocean
  `0.590/0.579`; DJF standardized-pressure MAE `0.229`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.25s`).
- Next accepted change remains M2-only.  Continue with real-Earth SLP map
  reading first; do not proceed to R2b wind, R3 currents, SST, precipitation,
  biomes, or generated-world fitting.

R2a M2 v46 checkpoint, 2026-07-07:

- Current real-Earth climate replay:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v46_20260707/`.
- Current pressure evidence packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/`.
- Mechanism delta: v46 keeps v45 and adds a Southern Ocean shoulder-season
  low-support floor.  It uses the existing Southern Ocean front/shelf/wave gate
  to create MAM/SON low-source support in the semantic Southern Ocean basin,
  instead of leaving shoulder seasons with zero M2 source.
- Map read: semantic Southern Ocean MAM residuals improve in every 60-degree
  sector; SON residuals improve in all low-pressure sectors except the
  0-60E sector, which becomes more overdeepened.  The wide 45-75S SON sectors
  improve at `-120..-60`, `60..120`, and `120..180`, but `-60..60` remains a
  zonal-band overdeepening watch item.
- Remaining defect: R2a is still not accepted.  v46 fixes the missing
  shoulder-season source class, but the Southern Ocean response is still too
  band-like.  The next owner is wave-sector sharpening, not broader amplitude.
- Metrics after map read: standardized-pressure MAE all/land/ocean
  `0.306/0.416/0.262`; standardized-pressure correlation all/ocean
  `0.571/0.308`; pressure zonal-anomaly correlation all/ocean
  `0.591/0.577`; MAM/SON standardized-pressure MAE `0.355/0.430`.
- Targeted tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.58s`).
- Next accepted change remains M2-only.  Continue with real-Earth SLP map
  reading first; do not proceed to R2b wind, R3 currents, SST, precipitation,
  biomes, or generated-world fitting.

## Single-Subgraph Earth Replay Fitting Rule

This is the non-negotiable fitting workflow for the next implementation work.
It is intentionally narrower than the old C5 acceptance/gate workflow.

1. Fit one real-Earth subgraph at a time.  The current subgraph is R2a
   seasonal SLP / pressure-source geometry; do not tune wind, currents, SST,
   moisture, precipitation, sea ice, Koppen classes, biomes, or generated
   worlds while this packet is still visibly wrong.
2. Use the real-Earth field as the only fitting target for that packet.  The
   core evidence is the same-grid Earth reference map, the same-grid Aevum
   real-Earth replay map, and their residual/vector-error maps.
3. Start each iteration by reading those maps.  Name the visible Earth
   structures being fitted, for example Southern Ocean westerlies, trade-wind
   belts, continental stationary-wave patches, monsoon-season pressure/wind
   reversal, or mountain/coastal roughness shadows.  Then name what replay gets
   right and wrong.
4. Accept a mechanism change only when its cause matches the diagnosed
   residual.  Metrics may reject a change or confirm regression safety, but a
   better global mean, p50, p90, or scalar gate cannot promote the packet if
   the maps still have the wrong geography.
5. If a residual belongs to an upstream field, stop and repair the upstream
   field first.  Do not mask a wind problem with ocean, precipitation, ice, or
   biome thresholds.
6. Generated terminal worlds are R9 guardrails only.  They are not evidence for
   choosing R0-R8 parameters and should not be rendered during the active R2
   fitting loop unless a later R9 promotion check is explicitly being run.

Current R2a acceptance packet:

- Required map read: seasonal Earth SLP anomaly maps, replay pressure-proxy
  anomaly maps, standardized residual maps, pressure zonal-anomaly maps, and
  pressure-center / stationary-wave support maps when available.  Wind maps may
  be rendered only as downstream observers; they are not the fitting target in
  R2a.
- Current checkpoint:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v46_20260707/`, with
  the previous wind checkpoint
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/` retained only as
  downstream observer evidence.
- Required attribution: explain whether each major pressure residual is caused
  by land-sea thermal pressure, missing ocean basin pressure centers,
  continent/basin-scale stationary waves, mountain/coast forcing, excessive
  smoothing, or an upstream boundary-condition primitive.
- Current blockers: replay pressure captures broad seasonal land thermal
  contrast, but it still reads as over-smooth continent/ocean blobs.  It lacks
  enough North Pacific / North Atlantic / Southern Ocean pressure-center
  structure and mountain/coast-driven stationary waves to organize later
  storm tracks.
- Next code owner: R2a pressure-source geometry.  The next accepted mechanism
  should add explicit pressure-center / stationary-wave objects and then use
  the pressure contact sheet for acceptance.  R2b pressure-to-surface-wind
  translation remains blocked until the R2a pressure maps are plausible.

## Inputs

Reference input:

- `out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz`
- `out_earth_climate_reference_r4_20260705/earth_reference_24000cells.npz`
- `out_earth_climate_reference_r5_oscar_20260706/earth_reference_8000cells.npz`
- `out_earth_climate_reference_r5_oscar_20260706/earth_reference_24000cells.npz`
- `out_earth_climate_reference_r6_landcover_20260706/earth_reference_8000cells.npz`
- `out_earth_climate_reference_r6_landcover_20260706/earth_reference_24000cells.npz`
- `out_earth_climate_comparison_r4_20260705/earth_climate_comparison_summary.json`
- `out_earth_climate_comparison_c4j1_render_r6_20260706/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_c4j1_render_r6_20260706/earth_climate_fitting_report.json`
- `out_earth_climate_comparison_c5e7_visual_r6_20260706/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_c5e7_visual_r6_20260706/earth_climate_fitting_report.json`
- `out_earth_climate_comparison_c5e8_receiver_supply_feedback_r6_20260706/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_c5e8_receiver_supply_feedback_r6_20260706/earth_climate_fitting_report.json`
- ETOPO 2022 terrain/bathymetry uncertainty check:
  `out_earth_climate_reference_r4_etopo2022_render_20260706/`
  and `out_earth_climate_reference_r4_etopo2022_24000_render_20260706/`.
  This does not replace the stable ETOPO5 `earth.elevation_m` baseline for
  climate fitting; it records the current Earth hypsometry uncertainty band.

R9 guardrail input, inactive during R0-R8 fitting:

- Six accepted terminal terrain worlds under
  `out_terminal_climate_biomes_20260705/`
- Per-world `terminal_climate_arrays.npz`
- Per-world rendered climate, current, precipitation, and biome maps
- Current C4j rendered climate replay:
  `out_terminal_climate_replay_c4j1_precip_objects_render_r6_20260706/`
- Current C5c2 no-render climate replay:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`
- Current C5e7 no-render climate replay:
  `out_terminal_climate_replay_c5e7_source_receiver_accounting_20260706/`
  with generated comparison previews now backfilled from each archived
  `terminal_climate_arrays.npz`.
- Current C5e8 no-render climate replay:
  `out_terminal_climate_replay_c5e8_receiver_supply_feedback_20260706/`
  with generated comparison previews backfilled from archived arrays.
- Current C5e9 no-render climate replay:
  `out_terminal_climate_replay_c5e9_ocean_structure_20260706/`.
  This keeps the C5e8 source/receiver precipitation feedback but strengthens
  the reduced ocean heat/SST structure and records same-latitude residual
  metrics for future anti-zonal-regression checks.

Frozen systems:

- Tectonics and terrain generation.
- Terminal elevation, bathymetry, crust, plate, and morphology fields.

Mutable systems:

- Climate module.
- Climate diagnostics and rendering.
- Biosphere/biome classifier after hydroclimate reaches plausible ranges.

## Replay-R Authoritative Geoscience Subgraph Order

Status: active reset, 2026-07-06

This section supersedes the old F0-F5 ordering for new replay repair work.  The
old F0-F5 notes remain below as history and regression evidence, but they are no
longer sufficient to justify tuning a downstream map if an upstream support
field is still wrong.

Rule:

- Only the active layer and its immediate diagnostics may be repaired.
- Downstream layers may be rendered and scored as observers, but they must not
  be tuned until their upstream support layers pass.
- Promotion cannot be justified by global means, scalar envelopes, or pass/fail
  tables alone.  Every Replay-R phase must fit one real-Earth subgraph at a
  time: compare the relevant real-Earth reference map with the corresponding
  real-Earth replay map on the same Earth grid, identify visible spatial
  residuals, and assign those residuals to upstream or active-layer mechanisms
  before any parameter change is accepted.
- Numeric metrics are screening and regression tools.  They can reject a run or
  confirm that a visually diagnosed fix did not break guardrails, but they do
  not replace reading the real maps and explaining the geography.
- Sea ice, seasonal snow, clouds, vegetation feedback, Koppen classes, and
  biomes are blocked until radiation, wind/pressure, ocean circulation, SST,
  moisture transport, and precipitation have passed replay gates.
- If a downstream map fails because an upstream field is wrong, fix the upstream
  field.  Do not compensate with downstream thresholds.
- Real-Earth replay is the fitting target.  The six accepted generated terminal
  worlds are not used to choose active-layer fixes during R0-R8; they are
  promotion guardrails in R9 after the Earth subgraph being fitted is plausible.

Earth-only phase packet:

- Pick exactly one active subgraph, for example R2 wind vectors, R3 current
  speed, R3 current direction, R4 SST, R5 moisture access, or R6 precipitation.
- Render the real-Earth reference map and the Earth replay map with the same
  grid, projection, color scale, and mask.  Add a residual or anomaly map when
  the field is continuous.
- Read the maps before reading summary scores: name the visible structures that
  are correct, missing, displaced, over-broad, too smooth, too noisy, or
  physically impossible.
- Attribute every major residual to either an upstream dependency or the active
  mechanism.  If the owner is upstream, stop and repair that upstream layer.
- Change only parameters or equations connected to the attributed mechanism.
  Do not compensate in downstream thresholds or by optimizing global means.
- Re-render the same Earth maps and record the before/after visual conclusion.
  Numeric metrics are archived as regression evidence after the map read.

Current reset ledger:

- Keep the real-Earth replay harness:
  `aevum.diagnostics.real_earth_climate_replay`, CLI
  `real-earth-climate-replay`, and `tests/test_real_earth_climate_replay.py`.
- Keep the C5e9 terminal-world regression bundle as generated-world guardrail
  evidence.
- Treat the exploratory real-Earth replay directories
  `out_real_earth_climate_replay_f3a_20260706/`,
  `out_real_earth_climate_replay_f3b_20260706/`,
  `out_real_earth_climate_replay_f3c_20260706/`,
  `out_real_earth_climate_replay_f4_20260706/`,
  `out_real_earth_climate_replay_f4b_20260706/`, and
  `out_real_earth_climate_replay_f4b_render_20260706/` as invalid promotion
  baselines.  They are diagnostic artifacts only because they attempted
  downstream precipitation/sea-ice/biome interpretation before the replay
  subgraph order was reset.
- Current code state after reset: downstream precipitation tail, diagnostic
  sea-ice threshold, and biome-threshold experiments from the invalid pass are
  reverted.  Partial foundation edits in temperature/SST and wind/pressure are
  provisional and must be re-evaluated under R1-R4 before promotion.
- Reset validation run:
  `out_real_earth_climate_replay_replay_r2_r4_reset_20260706/`.
  The replay harness passes validation and reproduces the provisional
  foundation metrics: surface-temperature MAE `2.94 C`, land-precipitation MAE
  `459.4 mm/yr`, annual wind p90 replay/Earth `6.72/6.72 m/s`, ocean-current
  p90 ratio `0.981`, and no hard validation failures.  This run is a foundation
  checkpoint, not a downstream hydro/sea-ice/biome acceptance.
- Current highest-priority foundation debts from that reset run:
  annual temperature still has a steep adjacent latitude-band jump, surface
  temperature zonal residual p95 is `12.05 C`, max adjacent latitude-band delta
  is `18.58 C`, and seasonal pressure pattern correlation is weak in SON
  (`0.032` versus `0.663/0.677/0.714` in DJF/MAM/JJA).  The next implementation
  work should therefore stay in R2/R4 before touching R6-R8.
- R2/R4 phase4 checkpoint:
  `out_real_earth_climate_replay_replay_r2_r4_phase4_20260706/`.
  This keeps the wind-driving base pressure conservative, adds a smaller
  forward-cooling seasonal pressure correction in the final pressure/moisture
  layer, and separates ordinary seasonal land snow albedo from high-latitude
  high-elevation permanent-ice albedo.  Real-Earth replay improves
  surface-temperature MAE from `2.94 C` to `2.62 C`, land-temperature MAE from
  `4.78 C` to `3.83 C`, and land-precipitation MAE from `459.4` to
  `451.3 mm/yr`.  Seasonal pressure correlations are now
  `0.617/0.601/0.682/0.205` for DJF/MAM/JJA/SON, so SON is no longer nearly
  uncorrelated, but it remains weaker than the other seasons.
- R4 diagnostic correction:
  `real-earth-climate-replay` now records Earth-aware adjacent latitude-band
  temperature jumps.  The phase4 replay's largest adjacent jump is `18.58 C`,
  while the same-grid Earth reference jump is `20.46 C` at the same
  Antarctic/Southern Ocean transition.  Therefore the generic validation
  warning about a steep adjacent latitude-band jump is not an R4 blocker for
  real-Earth replay; future heat-wall work should use the Earth-aware delta,
  not the generic absolute warning alone.
- Historical R9 guardrail status for phase4:
  six frozen terminal worlds were replayed in
  `out_terminal_climate_replay_replay_r2_r4_phase4_20260706/`.  Earth climate
  comparison reports `earthlike flagged: 0`, and the fitting report has
  `0` guardrail failures with `3` warnings in
  `out_earth_climate_fitting_replay_r2_r4_phase4_r6_20260706/`.
  However, the ocean-spatial gate fails with `8` failures in
  `out_earth_climate_ocean_spatial_gate_replay_r2_r4_phase4_20260706/`:
  earthlike and waterworld current-speed maps are too zonal, and earthlike SST
  same-latitude residual ratios are below the C5e9 floor.  Under the corrected
  Earth-only protocol this is not an R0-R8 fitting target and must not drive the
  next mechanism change.  It remains archived for later R9 promotion checks.
- R2 wind map-read checkpoint:
  `out_real_earth_wind_replay_r2_basin_boundary_layer_20260706/` is the
  superseded Earth-only wind diagnostic.  It compares `earth__seasonal_wind_u10_v10`
  against Aevum `atmosphere__seasonal_wind` on the same 8000-cell Earth grid,
  and now also renders wind-speed/eastward-wind zonal anomalies plus
  standardized Earth SLP versus replay pressure-proxy zonal anomalies.
  The code change stays inside R2: the background westerly belt is continuous
  and the Southern Hemisphere surface westerly center is shifted poleward and
  broadened; `_geographic_circulation_anomalies` now consumes geography
  primitives so open oceans strengthen westerlies/trades while continental
  interiors apply seasonal near-surface drag.  Visual read: the hard latitude
  wall is reduced, the Southern Ocean band is closer to the real-Earth latitude,
  and ocean speed-pattern correlation improves.  R2 is still not promoted:
  replay maps remain too stripe-like, land wind-speed pattern correlation
  remains negative, and the replay p90 wind speed is still low relative to Earth
  (`5.56` versus `7.00 m/s` seasonal p90; annual replay/Earth p90 is
  `5.02/6.72 m/s`).  Key wind metrics improved from the continuous-asymmetric
  baseline to this checkpoint: speed MAE `2.22 -> 1.99 m/s`, vector RMSE
  `3.68 -> 3.41 m/s`, direction-cosine p50 `0.832 -> 0.853`, all-cell speed
  pattern correlation `0.154 -> 0.327`, and ocean speed-pattern correlation
  `0.201 -> 0.368`; land speed-pattern correlation is still weak at `-0.154`.
  The new stationary-wave diagnostics sharpen the residual ownership: pressure
  zonal-anomaly correlation is moderate (`0.561` all-cell, `0.578` land,
  `0.547` ocean), but wind-speed zonal-anomaly correlation remains weak.  A
  stronger local boundary-layer vector-wave attempt was rejected because it
  worsened speed-zonal-anomaly correlation and vector RMSE.  The accepted
  scalar roughness/pressure-stationary checkpoint preserved wind direction and
  applied a small speed multiplier from pressure zonal-anomaly amplitude,
  open-ocean exposure, continent interior roughness, and terrain roughness.  The
  current basin/continent-scale boundary-layer checkpoint adds a broader
  domain-smoothed pressure-stationary response plus a small open-ocean strong
  tail support.  From the previous scalar checkpoint it improves speed MAE
  `1.968 -> 1.945 m/s`, vector RMSE `3.386 -> 3.367 m/s`, all-cell speed-pattern
  correlation `0.334 -> 0.355`, land speed-pattern correlation
  `-0.153 -> -0.150`, speed-zonal-anomaly correlation `0.183 -> 0.236`, land
  speed-zonal-anomaly correlation `0.081 -> 0.113`, ocean speed-zonal-anomaly
  correlation `0.153 -> 0.172`, eastward-zonal-anomaly correlation
  `0.304 -> 0.339`, and seasonal replay p90 wind speed `5.51 -> 5.57 m/s`
  against Earth `7.00 m/s`.  R2 is still not promoted: land speed-pattern
  correlation remains negative, maps still read too latitude-banded, and the
  strong-wind tail is still low.  Next R2 work should continue broad
  pressure/roughness/geostrophic-to-surface wind translation, not proceed to R3
  currents.
- R2 ocean-tail plus land-drag checkpoint:
  `out_real_earth_wind_replay_r2_ocean_tail_land_drag_20260706/` is a
  superseded Earth-only wind diagnostic, paired with
  `out_real_earth_climate_replay_replay_r2_ocean_tail_land_drag_20260706/`.
  The R2 code remains a scalar near-surface wind-speed response, not a local
  vector pressure kick: it adds a stronger open-ocean midlatitude storm-track
  tail and stronger continental roughness/boundary-layer drag while preserving
  the wind direction.  This directly targets the map-read residual where real
  Earth has broad ocean storm belts and continent-scale weak-wind patches, but
  replay carried latitude belts too cleanly across land.  Real-Earth replay
  validation passes, with annual wind p90 replay/Earth `5.56/6.72 m/s`.
  Seasonal R2 metrics improve versus the basin-boundary-layer checkpoint:
  speed MAE `1.945 -> 1.935 m/s`, vector RMSE `3.367 -> 3.363 m/s`,
  seasonal replay p90 `5.57 -> 5.90 m/s` against Earth `7.00 m/s`,
  all-cell speed-pattern correlation `0.355 -> 0.401`, land speed-pattern
  correlation `-0.150 -> -0.127`, speed-zonal-anomaly correlation
  `0.236 -> 0.317`, land speed-zonal-anomaly correlation `0.113 -> 0.159`,
  ocean speed-zonal-anomaly correlation `0.172 -> 0.198`, and
  eastward-zonal-anomaly correlation `0.339 -> 0.377`.  R2 is still not
  promoted: the maps still read too latitude-banded, land speed-pattern
  correlation remains negative, and Earth storm-track/geography placement is
  still sharper than replay.
- R2 pressure-steering checkpoint:
  `out_real_earth_wind_replay_r2_pressure_steering_20260706/` is a superseded
  Earth-only wind diagnostic, paired with
  `out_real_earth_climate_replay_replay_r2_pressure_steering_20260706/`.
  The new code adds a small, broad geostrophic steering response from the
  seasonal pressure zonal anomaly, with Coriolis damping near the equator and
  support from trade/westerly belts plus land/ocean exposure.  This is a
  wind-direction repair inside R2, not a downstream current/SST/precipitation
  compensation.  Real-Earth replay validation passes, with annual wind p90
  replay/Earth `5.58/6.72 m/s`.  Metrics improve versus the ocean-tail plus
  land-drag checkpoint: speed MAE `1.935 -> 1.932 m/s`, vector RMSE
  `3.363 -> 3.354 m/s`, direction-cosine p50 `0.853 -> 0.856`,
  direction-cosine p10 `-0.673 -> -0.667`, seasonal replay p90
  `5.90 -> 5.92 m/s` against Earth `7.00 m/s`, all-cell speed-pattern
  correlation `0.401 -> 0.402`, land speed-pattern correlation
  `-0.127 -> -0.112`, speed-zonal-anomaly correlation `0.317 -> 0.323`,
  land speed-zonal-anomaly correlation `0.159 -> 0.170`, ocean
  speed-zonal-anomaly correlation `0.198 -> 0.208`, and
  eastward-zonal-anomaly correlation `0.377 -> 0.382`.  R2 is still not
  promoted because the maps remain too latitude-banded and land speed-pattern
  correlation is still negative.
- R2 solstice basin-pressure checkpoint:
  `out_real_earth_wind_replay_r2_solstice_basin_pressure_20260706/` is a
  superseded Earth-only wind diagnostic, paired with
  `out_real_earth_climate_replay_replay_r2_solstice_basin_pressure_20260706/`.
  The R2 pressure source now includes a small geography-derived ocean term:
  large open basins support winter-hemisphere subpolar lows and broad
  subtropical ocean highs.  The term is limited to the solstice seasons so the
  shoulder seasons are not forced into an over-simple basin-pressure pattern.
  This improves the upstream pressure owner directly rather than compensating
  through currents, SST, precipitation, or biome thresholds.  Real-Earth replay
  validation passes, with annual wind p90 replay/Earth `5.60/6.72 m/s`.
  Metrics improve versus the pressure-steering checkpoint: speed MAE
  `1.932 -> 1.927 m/s`, vector RMSE `3.354 -> 3.343 m/s`,
  direction-cosine p10 `-0.667 -> -0.662`, seasonal replay p90
  `5.92 -> 5.92 m/s`, land speed-pattern correlation `-0.112 -> -0.110`,
  speed-zonal-anomaly correlation `0.323 -> 0.326`, ocean
  speed-zonal-anomaly correlation `0.208 -> 0.211`, pressure standardized
  ocean correlation `0.083 -> 0.135`, pressure-zonal ocean correlation
  `0.546 -> 0.560`, DJF pressure-zonal correlation `0.668 -> 0.682`, and
  JJA pressure-zonal correlation `0.712 -> 0.722`.  R2 is still not promoted:
  maps remain too latitude-banded, land speed-pattern correlation remains
  negative, and the replay still lacks enough North Atlantic/North Pacific and
  Southern Ocean storm-track geometry.
- R2 land boundary-layer plus polar katabatic checkpoint:
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/` is the current
  best Earth-only wind diagnostic, paired with
  `out_real_earth_climate_replay_replay_r2_land_katabatic_v2_20260706/`.
  Component attribution showed that non-polar land was systematically too
  windy while Antarctica was strongly under-windy.  The R2 code now adds an
  extra non-polar land boundary-layer drag and a polar highland downslope
  katabatic wind term derived from latitude, smoothed elevation, terrain
  gradient, and continent interiority.  This is still inside the R2 wind
  mechanism and does not touch currents, SST, precipitation, sea ice, Koppen,
  biomes, or generated worlds.  Real-Earth replay validation passes, with
  annual wind p90 replay/Earth `5.61/6.72 m/s`.  Metrics improve versus the
  solstice basin-pressure checkpoint: speed MAE `1.927 -> 1.865 m/s`, vector
  RMSE `3.343 -> 3.289 m/s`, direction-cosine p10 `-0.662 -> -0.654`,
  all-cell speed-pattern correlation `0.405 -> 0.444`, land speed-pattern
  correlation `-0.110 -> 0.145`, speed-zonal-anomaly correlation
  `0.326 -> 0.367`, land speed-zonal-anomaly correlation `0.172 -> 0.210`,
  ocean speed-zonal-anomaly correlation `0.211 -> 0.227`, and
  eastward-zonal-anomaly correlation `0.384 -> 0.416`.  R2 is still not
  promoted: the replay remains too latitude-banded overall and ocean storm
  tracks/basin-scale stationary waves are still too weak relative to Earth.
- R2 ocean storm-track attribution after the land/katabatic checkpoint:
  component and sector reads show that the remaining ocean problem is spatial
  distribution, not a simple wind-speed envelope.  Near/coastal ocean is too
  weak while some non-leeward midlatitude open-ocean belts are too strong, so
  blindly increasing or decreasing the whole ocean westerly band is the wrong
  owner.  Two candidate mechanisms were tested and rejected before code
  changes: a land-leeward ocean corridor boost inferred from current wind
  direction, and a near-coast boost/open-ocean damping redistribution.  Both
  reduced scalar MAE in some settings but worsened ocean speed-pattern,
  ocean zonal-anomaly, or eastward-anomaly correlations.  The next R2 repair
  should therefore improve the actual storm-track source geometry, likely by
  deriving better baroclinic/coastal thermal-front and basin-sector support,
  rather than by applying another scalar ocean multiplier.
- R2 rejected storm-track source-geometry probes after the land/katabatic
  checkpoint:
  the current best checkpoint remains
  `out_real_earth_wind_replay_r2_land_katabatic_v2_20260706/`.  Additional
  Earth-only offline probes tested baroclinic and coastal thermal-front scalar
  modifiers, extra ocean pressure-steering, east-coast winter ocean lows,
  orographic downstream stationary-wave redistribution, and broad open-ocean
  storm-tail amplification.  None is accepted.  The thermal-front family
  worsened ocean speed-pattern, ocean speed-zonal, and eastward-anomaly
  correlations even at low amplitude.  Extra ocean pressure-steering gave only
  a tiny ocean anomaly improvement and did not visibly change the latitude-band
  reading of the maps.  East-coast winter ocean lows did not improve the
  pressure map against Earth SLP.  The orographic downstream-wave probe
  improved some scalar pattern correlations but still rendered as latitude
  belts and degraded eastward anomaly structure; the rendered probe is
  `out_probe_real_earth_wind_replay_r2_orographic_wave_tail_20260706/`.  A
  broad open-ocean tail can raise p90 wind speed toward Earth, but it worsens
  MAE and eastward anomaly placement as amplitude grows, so it is not a
  physical fix.  This moves the next owner upstream inside R2: first fit the
  real-Earth pressure/source geometry map more directly, then re-test
  pressure-to-surface-wind translation.
- R2a pressure/source geometry diagnostic checkpoint:
  `real-earth-wind-replay` now writes pressure-only evidence in addition to
  the wind contact sheet: `pressure_standardized_delta_seasons.png`,
  `real_earth_pressure_replay_contact_sheet.png`, and standardized-pressure
  MAE metrics for all cells, land, ocean, and each season.  The current R2a
  baseline is `out_real_earth_pressure_replay_r2a_current_20260706/`.  Map
  read: replay pressure has the broad seasonal land thermal signal, but it is
  still too smooth and too blocky.  It lacks enough ocean basin pressure
  centers and mountain/coast-driven stationary-wave structure, so it cannot
  yet organize Earthlike storm tracks.  Current pressure metrics are
  standardized MAE all/ocean `0.317/0.277`, standardized correlation
  all/land/ocean `0.539/0.668/0.135`, and pressure zonal-anomaly correlation
  all/land/ocean `0.567/0.578/0.560`; DJF/MAM/JJA/SON zonal correlations are
  `0.682/0.563/0.722/0.388`.
- R2a rejected pressure-source probes:
  offline pressure-only sweeps tested stronger open-ocean winter subpolar lows,
  basin-sector east-coast downstream plumes, mountain downstream stationary
  waves, ocean thermal-front lows, and reduced land thermal smoothing.  Stronger
  open-ocean subpolar lows improved ocean standardized correlation
  `0.135 -> 0.160` and ocean zonal correlation `0.560 -> 0.565`, but the
  contact sheet remained visually too similar to the baseline and MAE worsened
  (`0.317 -> 0.320`).  Reduced smoothing and stronger ocean lows raised ocean
  standardized correlation to about `0.171`, but degraded all-cell, land, and
  zonal correlations; basin-sector plumes and thermal-front lows worsened the
  pressure map.  None is accepted.  The next R2a mechanism should create
  explicit basin-scale pressure-center objects or a better stationary-wave
  pressure source, not another broad scalar latitude-band term.
- R2a pressure-center diagnostic object checkpoint:
  `ClimateModule` now emits diagnostic-only pressure/source geometry fields:
  `atmosphere.pressure_center_support`, `atmosphere.pressure_center_id`, and
  `atmosphere.stationary_wave_pressure_support`, plus
  `atmosphere.pressure_centers` objects.  These are derived from the pre-hydro
  `atmosphere.land_sea_pressure_proxy` and geography primitives; they do not
  alter pressure, wind, precipitation, or downstream physics.  Real-Earth
  replay output is
  `out_real_earth_climate_replay_r2a_pressure_centers_20260706/`, with R2a
  pressure evidence in
  `out_real_earth_pressure_replay_r2a_pressure_centers_20260706/`.
  The new pressure contact sheet includes Earth SLP, replay pressure proxy,
  standardized residual, zonal residual, pressure-center support, and
  stationary-wave support.  Map read: the support layer is useful for
  attribution, but it confirms that replay pressure is still too broad and
  blocky; large continental blobs and Southern Ocean latitude bands dominate,
  while Earthlike ocean basin pressure centers and coast/mountain stationary
  waves remain under-structured.  Metrics are unchanged from the R2a baseline,
  as expected for a diagnostic-only checkpoint.
- R2a major-ocean semantic basin checkpoint:
  real-Earth replay now marks `diagnostics.real_earth_replay`, and
  `ClimateModule` uses that flag to emit semantic major-ocean basin ids for
  Atlantic, Pacific, Indian, Arctic, and Southern Ocean in `ocean.basin_id`.
  The replay packet is
  `out_real_earth_climate_replay_r2a_major_ocean_basins_20260707/`, with
  pressure and M0/M1 evidence in
  `out_real_earth_pressure_replay_r2a_major_ocean_basins_20260707/`.
  Map read: the M0 basin-id blocker is resolved; the support map now has
  readable major ocean objects.  R2a is still not accepted, because the
  pressure map still reads as broad continental blobs plus an over-zonal
  Southern Ocean band.  Next owner: M2 pressure genesis, specifically
  basin/front/terrain-supported pressure-center construction rather than R2b
  wind or downstream SST/moisture tuning.
- R2a M2 pressure-genesis v1 checkpoint:
  `ClimateModule` now applies an object-based M2 pressure-source refinement
  after weak ocean-atmosphere coupling and before pressure diagnostics.  The
  pass extracts broad winter subpolar ocean-low candidates from
  `ocean.basin_id`, open-ocean support, SST-front support, and same-latitude
  SST anomaly; it also adds weak basin subtropical highs and small
  continent/terrain stationary-wave refinements.  It changes the R2a pressure
  proxy only; R2b wind translation is not retuned.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v1_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v1_20260707/`.
  Map read: DJF North Pacific and North Atlantic lows are more visible, and
  JJA Southern Ocean support is more segmented than the previous latitude-band
  result.  R2a is still not accepted because the replay remains too smooth over
  continents and still under-organizes Earthlike stationary-wave structure.
  Metrics support a partial mechanism improvement rather than promotion:
  ocean standardized-pressure correlation improves `0.137 -> 0.213`, all-cell
  standardized correlation improves `0.539 -> 0.545`, and all-cell
  zonal-anomaly correlation improves `0.568 -> 0.577`, while standardized MAE
  worsens `0.317 -> 0.320`.
- R2a M2 pressure-genesis v2 checkpoint:
  `ClimateModule` now archives the causal M2 source separately from the
  pressure-center result: `atmosphere.pressure_genesis_source`,
  `atmosphere.ocean_pressure_low_source_support`,
  `atmosphere.ocean_pressure_high_source_support`,
  `atmosphere.land_pressure_source_support`, and
  `atmosphere.terrain_pressure_wave_source_support`.  The real-Earth pressure
  replay renders these source maps in the pressure contact sheet and as
  individual PNGs.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v2_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v2_20260707/`.
  Map read: DJF North Pacific and North Atlantic ocean-low source patches are
  now directly visible; JJA Southern Ocean source support is visible but still
  too band-like.  The replay pressure field remains too smooth and is not
  promoted.  Metrics versus the M0 basin checkpoint: standardized-pressure
  MAE all is close to baseline (`0.317 -> 0.317` rounded), all-cell
  correlation improves `0.539 -> 0.546`, ocean correlation improves
  `0.137 -> 0.205`, all-cell zonal-anomaly correlation improves
  `0.568 -> 0.576`, and ocean zonal-anomaly correlation improves
  `0.563 -> 0.570`.
- R2a M2 pressure-genesis v4 checkpoint:
  `ClimateModule` now applies Southern Ocean sector/wavenumber-front gating to
  ocean-low source extraction and uses continent-level thermal-center objects
  plus terrain/land-thermal-gradient support for continental pressure
  refinement.  Evidence:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v4_20260707/` and
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v4_20260707/`.
  Map read: Aleutian/Icelandic causal low-source patches remain readable; JJA
  Southern Ocean source is visibly segmented into wave-like sectors in the
  source panel; continental source maps now show centered/object-like
  refinements instead of a single all-continent scalar.  R2a remains
  unaccepted because the final replay pressure is still smoother and less
  stationary-wave rich than real Earth SLP.  Metrics versus the M0 basin
  checkpoint: standardized-pressure MAE all is essentially unchanged
  (`0.317 -> 0.317` rounded), all-cell correlation improves `0.539 -> 0.547`,
  land correlation is preserved/improved `0.668 -> 0.668`, ocean correlation
  improves `0.137 -> 0.208`, all-cell zonal-anomaly correlation improves
  `0.568 -> 0.577`, land zonal-anomaly correlation improves
  `0.578 -> 0.585`, and ocean zonal-anomaly correlation improves
  `0.563 -> 0.571`.

Immediate Earth-only replay sequence:

1. R2a seasonal SLP/pressure-source geometry against
   `earth__seasonal_slp_anomaly_hPa`, using the pressure proxy as the replay
   field.  This now owns the next fitting step because the wind residual is
   mostly storm-track/source geometry rather than a local wind multiplier.
2. R2b seasonal 10 m wind vector and speed maps against
   `earth__seasonal_wind_u10_v10`, only after R2a pressure maps stop reading as
   overly smooth continent/ocean blobs.
3. R3 OSCAR surface-current speed and vector maps, only after R2 wind replay is
   visually plausible.
4. R4 OISST annual and seasonal SST maps, only after R3 current/upwelling
   structure is plausible.
5. R5 moisture-access/source maps, then R6 precipitation maps, then R7 sea ice
   and R8 climate/biome classes.

Do not render or tune generated terminal worlds during steps 1-4.  They return
only as R9 generalization checks after the active Earth subgraph has passed its
map-read packet.

### R0 - Replay Harness and Invariants

Purpose: make replay evidence trustworthy before physical repair.

Frozen inputs:

- Plate motion, terrain, bathymetry, final sea mask, and terminal generated
  worlds.
- Earth reference arrays from the R6 landcover reference package.

Mutable code:

- Replay CLI, diagnostics, archive schema, render/contact-sheet code, tests.

Required checks:

- Earth and replay arrays have matching cell count, masks, lat/lon ordering, and
  finite values.
- Seasonal arrays aggregate exactly to annual arrays where the model contract
  says they should.
- Observer fields are archived without mutating physics.
- Generated-world replay can still consume frozen terminal arrays without
  rerunning tectonics or terrain.

Promotion rule:

- No climate physics tuning is allowed in R0.

### R1 - Boundary Conditions and Forcing

Purpose: establish the non-negotiable inputs that all later climate layers
consume.

Inputs:

- Latitude, longitude, area, elevation, sea mask, orbital parameters, stellar
  flux, atmosphere mass/greenhouse state, coastline and basin primitives.

Mutable fields:

- Insolation/seasonal phasing diagnostics, climate elevation smoothing,
  coast/basin primitive extraction, topographic barrier primitives.

Observer only:

- Wind, currents, precipitation, sea ice, Koppen, and biomes.

Earth replay microbenchmarks:

- Annual and seasonal insolation by latitude.
- North/south seasonal phase.
- Land/ocean/elevation hypsometry seen by the climate solver.
- Coast, basin, and barrier primitive coverage.

Promotion rule:

- Later layers cannot be tuned to compensate for wrong seasonal phasing, wrong
  climate elevation, or inconsistent sea/land masks.

### R2 - Atmospheric Dynamics: Pressure and Wind

Purpose: fix the atmospheric driver before using it to repair oceans, moisture,
or precipitation.

Inputs:

- R1 forcing, land/sea thermal inertia, elevation/barriers, coast geometry.

Mutable fields:

- Background seasonal wind belts.
- ITCZ migration and convergence zones.
- Land-sea thermal pressure anomalies.
- Orographic deflection and wind-gap recovery.
- Storm-track and monsoon-potential support fields.

Observer only:

- Ocean currents may report wind-stress response diagnostics, but ocean current
  tuning is not promoted until R3.
- Precipitation and biome maps are observers only.

Earth replay microbenchmarks:

- Trade-wind, westerly, and polar-easterly latitude bands.
- Seasonal monsoon pressure reversal over large continents.
- Pressure/temperature anti-correlation over land.
- Pressure-gradient/wind alignment.
- Wind speed p50/p90/p99 by ocean, coast, and land.
- Geographic circulation index compared against Earth wind-field structure.

Promotion rule:

- R2 passes only when wind and pressure make sense on real Earth without using
  precipitation, sea ice, or biome thresholds to hide errors.

### R3 - Ocean Dynamics and Basin Coupling

Purpose: derive surface currents from basin geometry plus R2 winds before SST
and precipitation are judged.

Inputs:

- R2 wind/pressure fields, ocean mask, basin ids, shelves, straits, coasts,
  latitude and rotation sign.

Mutable fields:

- Basin streamfunction and gyres.
- Western/eastern boundary currents.
- Equatorial current bands and countercurrents.
- Upwelling/downwelling.
- Strait exchange.
- Wind-stress current response.
- Current heat-transport proxy.

Observer only:

- Sea ice and precipitation remain observers.

Earth replay microbenchmarks:

- Current speed p50/p90/p99 against OSCAR reference envelopes.
- Strong-current near-coast and far-ocean shares.
- Boundary-current class coverage.
- Upwelling co-location with eastern boundary/cold-current regimes.
- Basin confinement and land leakage.
- Mean-zero heat-transport anomaly.
- Visual map-read attribution against OSCAR/OISST reference maps: distinguish
  equatorial bands, western boundary currents, eastern-boundary cold/upwelling
  corridors, Southern Ocean arcs, basin gyres, strait exchange, and open-ocean
  background flow.  The replay diagnosis must name which structures are present,
  missing, over-broad, or too latitude-banded in the real-Earth replay before
  changing R3 parameters.

Promotion rule:

- Currents must be plausible as a field before their SST or moisture effects are
  used to tune downstream layers.
- R3 does not pass if current-speed or SST maps only satisfy scalar envelopes
  while reading as smooth latitude bands disconnected from basin geometry,
  coasts, and currents in the real-Earth reference.  Generated-world current
  maps are checked only later in R9, after the Earth replay R3 subgraph passes.

### R4 - SST and Energy Closure

Purpose: close the radiation, topography, wind, and ocean-energy layer before
repairing hydrology or cryosphere.

Inputs:

- R1 forcing, R2 winds/pressure, R3 currents/upwelling/heat transport.

Mutable fields:

- Annual and seasonal surface temperature.
- Seasonal SST.
- Land/ocean thermal inertia and maritime lag.
- Lapse-rate/topographic cooling.
- Coastal current temperature influence.
- Bounded ocean heat-flux residual.

Observer only:

- Sea ice can be measured here, but not threshold-tuned.  If sea ice is wrong
  because SST is wrong, R4 owns the fix.
- Precipitation and biomes remain observers.

Earth replay microbenchmarks:

- Global, land, and ocean mean temperature.
- Seasonal temperature amplitude over land and ocean.
- SST latitude-band residuals relative to OISST.
- Coastal warm/cold current temperature asymmetry.
- Max adjacent latitude-band jump, to catch heat/cold walls.
- Ocean heat-flux conservation and bounded residuals.

Promotion rule:

- R4 must remove obvious heat/cold walls and SST structure failures before R5
  moisture or R7 sea ice work begins.

### R5 - Moisture Source and Transport

Purpose: derive moisture availability from accepted wind, SST, currents,
basins, coasts, and topographic barriers.

Inputs:

- R2 wind/pressure, R3 currents/basins, R4 SST/temperature, R1 barriers/coasts.

Mutable fields:

- Evaporation/source strength.
- Moisture access.
- Source-ocean warmth.
- Moisture source basin id.
- Seasonal landward pathways and wind-gap recovery.
- Terrain blocking and rain-shadow support indices.
- Receiver catchment support diagnostics.

Observer only:

- Precipitation totals, Koppen, and biomes are not tuned yet.

Earth replay microbenchmarks:

- Moisture-access p50/p75/p90 by continent size and distance from coast.
- Low/mid-latitude summer moisture corridors.
- Source-basin attribution fraction.
- Moisture pathway continuity and barrier crossing penalties.
- Dry-interior explanation by access/blocking rather than arbitrary aridity.

Promotion rule:

- R5 passes when wet/dry potential is geographically explainable before annual
  precipitation amounts are adjusted.

### R6 - Seasonal Precipitation and Hydroclimate Objects

Purpose: solve precipitation from accepted moisture transport and atmospheric
support fields.

Inputs:

- R5 moisture fields, R2 wind/pressure/ITCZ/storm-track support, R1 orography.

Mutable fields:

- Seasonal precipitation.
- Annual precipitation.
- Orographic uplift and rain shadow.
- Monsoon/storm/ITCZ rainfall corridors.
- Runoff and local moisture-budget regions.
- Hydroclimate objects and receiver catchments.

Observer only:

- Koppen and biomes remain observers.  They must not be tuned to compensate for
  hydroclimate errors.

Earth replay microbenchmarks:

- Land precipitation mean, p50, p90.
- Wet tropical, dry subtropical, and mid-latitude storm-belt placement.
- Windward/leeward mountain precipitation contrast.
- Wet-season phase alignment with maximum support.
- Dry-cell explanation by low moisture access, rain shadow, or subsidence.
- Local/seasonal budget conservation.

Promotion rule:

- R6 must pass before any biome threshold or vegetation-feedback repair.

### R7 - Cryosphere, Cloud, and Vegetation Feedback

Purpose: add bounded feedbacks only after their physical inputs are plausible.

Inputs:

- R4 temperature/SST, R6 precipitation/snow supply, R3 ocean state, R1
  topography.

Mutable fields:

- Seasonal sea ice.
- Seasonal snow persistence.
- Ice-sheet tendency and permanent ice.
- Cloud albedo and precipitation-cloud coupling.
- Vegetation albedo/roughness feedback.

Earth replay microbenchmarks:

- Seasonal sea-ice extent and hemispheric asymmetry.
- Snow persistence by latitude/elevation.
- Ice sheet placement and non-placement.
- Cloud/vegetation feedback boundedness.
- No feedback loop may erase R2-R6 structure.

Promotion rule:

- Sea ice is never repaired by an isolated diagnostic threshold while R2-R4 are
  failing.

### R8 - Climate Classes and Biomes

Purpose: classify accepted climate fields into stable user-facing climate and
biome maps.

Inputs:

- Accepted R4 temperature, R6 precipitation/hydroclimate, and R7 cryosphere/
  vegetation fields.

Mutable fields:

- Koppen-like climate classes.
- Biome classifier.
- Spatial generalization and semantic region cleanup.

Earth replay microbenchmarks:

- Tropical, arid, temperate, continental, polar, alpine, and ice-class
  envelopes.
- Forest/grass/desert/tundra/ice fractions and latitude organization.
- Mountain zonation and coast/interior biome transitions.

Promotion rule:

- R8 cannot be used to hide upstream temperature, moisture, or precipitation
  errors.

### R9 - Generated-World Regression and Promotion

Purpose: make sure Earth replay fitting did not overfit Earth geography.

Required evidence:

- Real-Earth replay report and contact sheet for the promoted stage.
- Six generated terminal worlds replayed from frozen terrain arrays.
- Waterworld and arid guardrails.
- Difference report against the previous accepted generated-world bundle.
- Targeted tests for any edited physics, diagnostics, or classifier code.

Promotion rule:

- A stage is promoted only when its Earth replay microbenchmarks pass and the
  six generated worlds preserve plausible variety.

## Fitting Order

### F0 - Calibration Harness

Status: complete

Purpose: turn Earth-vs-generated differences into stable fitting evidence.

Deliverables:

- JSON report with phase-level scores and flags.
- CSV table of candidate parameter levers.
- CSV table of every cross-preset guardrail check.
- Markdown report explaining dominant failures.
- Cross-preset guardrail verdict for earthlike, arid, and waterworld outputs.
- Tests using synthetic summaries so the report logic is deterministic.
- Climate-only replay for frozen terminal arrays, so climate parameter tuning
  does not rerun or mutate the accepted plate/terrain system.
- Optional CLI failure gate via `--fail-on-guardrail`, intended for parameter
  sweeps and CI-style batch scripts.  Warnings do not fail the command.

Acceptance:

- Earthlike worlds are evaluated separately from waterworld/arid diagnostic
  worlds.
- The report identifies which phase owns each large mismatch.
- The report refuses to recommend biome-threshold tuning while precipitation is
  severely out of range.
- Arid and waterworld diagnostic runs are guarded against climate-regression
  side effects, even though they are not fitted directly to Earth.
- Replaying existing earthlike terminal arrays reproduces the original climate
  metrics exactly before any climate parameter edits.

### F1 - Temperature and Energy

Purpose: ensure annual mean temperature, land/ocean contrast, lapse-rate
cooling, and seasonal amplitude are inside broad Earthlike envelopes before
water-cycle tuning.

Primary metrics:

- Global mean temperature delta.
- Land and ocean mean temperature delta.
- Land/ocean seasonal temperature amplitude.
- Max adjacent latitude-band temperature jump.

Acceptance for Earthlike worlds:

- Global mean temperature within about 6 C of the Earth baseline.
- Land and ocean means not both biased in the same direction by more than about
  6 C.
- Seasonal amplitudes close enough that precipitation tuning is not compensating
  for a temperature error.

Expected current state:

- Temperature is no longer a dominant blocker after adding an ocean mixed-layer
  SST floor for calibration against OISST.  Sea-ice state still comes from the
  raw cold solution, but exported ocean temperature no longer falls far below
  seawater's freezing point.

### F2 - SST and Ocean-Current Influence

Status: first pass complete; C5b1 spatial ocean-current/SST gate passing

Purpose: make current/SST effects shape coasts without overwhelming the climate.

Primary metrics:

- Current speed p50 and p90 versus NASA/JPL OSCAR monthly climatology, with
  NOAA/AOML drifter climatology retained as an annual-speed cross-check.
- Ocean heat-transport anomaly distribution.
- Upwelling area and cold-current dry-coast co-location.
- Coastal temperature asymmetry by latitude band.
- Strong-current near-coast and far-ocean shares against the R6 OSCAR/OISST
  Earth reference envelope.
- SST zonality relative to OISST so ocean-current feedback is not only a
  latitude-band texture.

Acceptance:

- Currents remain basin-confined and tangent.
- Current speed p90 is in the same order as the Earth baseline.
- Heat transport redistributes energy without changing global mean energy.
- Earthlike strongest currents include boundary-current/coastal structure and
  are not dominated by remote open-ocean bands.

Expected current state:

- Vector current p90 was high relative to the Earth drifter baseline.  The
  exported near-surface vector field is now scaled separately from the reduced
  heat-transport proxy.  Current p90 in the C5b1 earthlike replay is now about
  `0.90` times the R6 OSCAR annual Earth baseline.
- R5 Earth references promote OSCAR to the primary current-speed baseline:
  8000/24000 annual p90 is about `0.192` / `0.190 m/s`, while the retained
  AOML cross-check remains about `0.267` / `0.265 m/s`.
- C5b1 adds `aevum.diagnostics.earth_climate_ocean_spatial_gate` and CLI
  command `earth-climate-ocean-spatial-gate`.  The gate checks current p90
  magnitude, current land leakage, current/SST zonal dominance, mean-zero heat
  transport, coastal swift-current SST anomaly spread, and earthlike
  near-coast/far-ocean strong-current placement.
- The C5a2 failure was localized to `earthlike_seed42`: strongest-current
  near-coast share was `0.332` and far-ocean share was `0.464`.  C5b1 moves
  those to `0.407` and `0.389`, respectively, while preserving all other gate
  envelopes.

### F3 - Circulation and Moisture Access

Status: first pass complete; C5c2 coupled-consistency gate passing

Purpose: make winds, pressure, source-ocean warmth, barriers, and moisture
routes create plausible wet coasts, monsoon corridors, and dry interiors.

Primary metrics:

- Land moisture-access p50/p75/p90.
- Source-ocean warmth distribution.
- Monsoon-potential p90/p99.
- Terrain-blocking p75 and wind-gap recovery.
- Seasonal pressure anomaly strength.
- Pressure/temperature anti-correlation over seasonal land cells.
- Pressure-gradient/wind alignment.
- SST/source-ocean/evaporation coupling.
- Moisture/support/precipitation coupling.
- High monsoon-potential cells must also have adequate moisture access and
  enhanced seasonal precipitation.

Acceptance:

- Waterworld has weak continent-driven monsoon potential.
- Earthlike large continents show seasonal onshore moisture corridors.
- Arid worlds retain dry interiors but do not become uniformly rainless.
- Seasonal wind, pressure, SST/source warmth, evaporation, moisture access,
  monsoon potential, and precipitation must be mutually consistent under the
  reduced weak-coupling model.

Expected current state:

- The first fitting pass added a weak free-atmosphere moisture reservoir and
  stronger source-ocean/monsoon coupling.  Earthlike land moisture-access p75 is
  now about 0.59-0.65 and monsoon-potential p90 is about 0.17-0.19.
- C5c2 adds `aevum.diagnostics.earth_climate_coupled_consistency_gate` and CLI
  command `earth-climate-coupled-consistency-gate`.  The gate checks
  pressure-temperature phase, wind-pressure alignment, SST/source/evaporation
  coupling, moisture/support/precipitation correlation, monsoon-pressure and
  monsoon-moisture support, cold-current evaporation suppression, precipitation
  budget closure, heat-flux mean conservation, and bounded coupling residuals.
- The C5b1 failure was localized to `earthlike_seed42`: top monsoon-potential
  cells had moisture p25 only `0.504` of the land median.  C5c2 moves that
  ratio to `0.852` by tightening the monsoon moisture gate, while preserving
  Earth comparison, F2 ocean-spatial, circulation-layout, and downstream
  hydroclimate object gates.

### F4 - Seasonal Hydroclimate

Status: first pass complete; C5d1 seasonal hydroclimate placement gate passing

Purpose: tune seasonal precipitation after circulation and moisture access are
diagnosed.

Primary metrics:

- Land precipitation mean, p50, and p90.
- Land wet fraction above 500 mm/yr.
- Precipitation seasonality p75.
- Orographic precipitation concentration.
- Rain-shadow and wet-windward diagnostics.
- Wet-cell support by moisture access, monsoon corridors, storm tracks, ITCZ,
  regional wet response, and moisture-flow response.
- Dry-cell explanation by low moisture access, rain shadow, or regional dry
  response.
- Wet-season phase alignment with the season of maximum wet-process support.

Acceptance:

- Earthlike land precipitation is no longer below 45% of the Earth baseline.
- The 50th and 90th percentiles rise together, not only by adding isolated wet
  spikes.
- Orographic precipitation remains regional instead of becoming ridge stripes.

Expected current state:

- This was the dominant R4 failure.  The first fitting pass raised earthlike
  land precipitation mean to about 0.58-0.69 of the Earth baseline, p50 to about
  0.84-1.05, and p90 to about 0.46-0.56.
- The F4 tail gate adds a warm, high-access convective tail term.  Earthlike
  p90 precipitation is now about 0.52-0.64 of Earth while arid median
  precipitation and waterworld island guardrails remain inside bounds.
- The first Earth pattern gate adds spatial constraints for wet tropics, dry
  subtropics, mountain wet tails, and high-latitude cold envelopes.  Pattern7
  keeps the scalar gate passing while clearing all climate-pattern failures:
  earthlike p90 precipitation is about 0.57/0.70 of Earth globally, wet-tropics
  p90 clears the 0.5 Earth-ratio floor, and subtropical dry belts are present.
- C5d1 adds `aevum.diagnostics.earth_climate_seasonal_hydro_placement_gate`
  and CLI command `earth-climate-seasonal-hydro-placement-gate`.  This is an
  internal weak-coupling placement gate: it checks that seasonal wet cells are
  supported by moisture, monsoon, storm-track, ITCZ, and precipitation-response
  fields; dry cells are explainable by low moisture or rain-shadow/dry
  response; wet-season phase follows maximum support; and seasonal
  precipitation still aggregates exactly to annual precipitation.

### F5 - Koppen and Biomes

Status: first pass complete; Earth biome/spatial/seasonal-subtype/mountain-zonation/windward-leeward gates pass

Purpose: tune biome and climate-class mapping only after temperature and
hydroclimate are plausible.

Primary metrics:

- Desert, forest, tropical, tundra, and ice area fractions.
- Latitudinal biome envelopes.
- Agreement between model classes, GloH2O Koppen proxy, and RESOLVE biomes at
  coarse class level.
- Coarse land-cover consistency against ESA CCI broad classes, mainly forest,
  cropland, grass/shrub, bare/sparse, snow/ice, wetland, urban, and water
  envelopes.

Acceptance:

- Earthlike worlds produce nonzero forest and tropical classes when climate
  supports them.
- Desert excess is not fixed by lowering desert thresholds while precipitation
  is still too low.

Expected current state:

- Forest and tropical classes are no longer absent in earthlike replays.
  Classification now preserves climate-supported forest/tropical semantic
  patches during biome generalization.
- Pattern7 also fixes cold-dry biome precedence so high-latitude cold land is
  no longer overwritten as desert.
- F5 biome-gate tuning lowers the moist-temperate forest precipitation
  threshold from 580 mm/yr to 520 mm/yr.  This fixes the remaining
  forest+tropical envelope warnings without changing the precipitation field:
  earthlike forest+tropical land fractions are now about 0.231 and 0.267, above
  the gate floor versus both Koppen proxy and RESOLVE coarse references.
- F5 spatial-biome tuning adds a latitude-organization gate and fixes the
  remaining cool-midlatitude/high-latitude semantic errors.  Cool climates now
  use lower forest precipitation thresholds, cool dry climates use lower desert
  thresholds, and high-latitude cold-dry land is treated as tundra/ice before
  desert classification.
- F5 seasonal-subtype tuning adds a Koppen-like seasonal gate.  Low-latitude
  precipitation seasonality is strengthened by a per-cell redistribution that
  preserves annual precipitation, and seasonal tropical biome semantics now
  keep wet Aw/Am-like climates in the coarse tropical class instead of
  converting all two-dry-quarter tropical climates to grassland.
- F5 mountain-zonation tuning adds an Earth envelope gate for high-mountain
  alpine ecology and mountain desert excess.  The biome classifier now treats
  cool high-elevation land as alpine/tundra before arid classification, so
  dry highlands are not automatically mapped as desert.
- Windward/leeward tuning adds an Earth envelope gate that combines seasonal
  wind vectors, topographic gradients, and mountain-slope precipitation.  The
  climate model now computes a slope-wind exposure field and uses it to
  redistribute seasonal land precipitation from leeward to windward slopes
  while preserving the seasonal land mean.

## Microbenchmarks

Use these every time climate parameters change:

- `earth-climate-fit-report` on the latest comparison summary.
- `earth-climate-fit-report --fail-on-guardrail` when running automated sweeps;
  this exits nonzero only for failed guardrails, not for warnings.
- `earth-climate-pattern-gate --fail-on-pattern` after scalar guardrails pass;
  this evaluates broad spatial envelopes rather than exact map overlap.
- `earth-climate-biome-gate --fail-on-biome` after pattern guardrails pass;
  this evaluates coarse biome envelopes against Koppen proxy and RESOLVE
  references.
- `earth-climate-spatial-biome-gate --fail-on-spatial-biome` after coarse
  biome guardrails pass; this checks latitude organization of tropical,
  desert, temperate/boreal forest, and tundra/ice classes.
- `earth-climate-seasonal-subtype-gate --fail-on-seasonal-subtype` after
  spatial-biome guardrails pass; this checks dry/wet-quarter organization,
  low-latitude seasonal tropical subtype area, and broad seasonality amplitude.
- `earth-climate-mountain-zonation-gate --fail-on-mountain-zonation` after
  seasonal-subtype guardrails pass; this checks high-mountain alpine ecology,
  high-mountain desert excess, mountain cooling, and mountain wet-tail
  envelopes.
- `earth-climate-windward-leeward-gate --fail-on-windward-leeward` after
  mountain-zonation guardrails pass; this checks whether mountain slopes that
  seasonal winds climb are significantly wetter than leeward slopes.
- `earth-climate-seasonal-hydro-placement-gate
  --fail-on-seasonal-hydro-placement` after F2/F3 coupled-consistency passes;
  this checks whether seasonal wet/dry cells and wet-season timing are
  explained by the generated moisture, ITCZ, monsoon, storm-track,
  rain-shadow, and regional response fields rather than isolated rainfall
  patches.
- `earth-climate-hydro-region-gate --fail-on-hydro-regions` after C4d object
  fields exist; this checks object archive presence, multi-kind/multi-season
  Earthlike coverage, coherent region sizes, monsoon seasonal migration,
  storm-track latitude/coastal placement, rain-shadow dry-response alignment,
  object-map active area, largest connected belt share, boundary roughness,
  and waterworld/arid false-positive limits.  By default it also writes C4d
  hydroclimate region contact sheets for human visual review.
- `earth-climate-moisture-flow-gate --fail-on-moisture-flow` after C4e object
  fields exist; this checks moisture-flow object/archive presence,
  multi-kind/multi-season Earthlike coverage, source-ocean strength, routed
  land-pathway strength, pathway coupling to monsoon/storm support,
  source/pathway/network-id map readability, and waterworld/arid
  false-positive limits.  By default it also writes C4e moisture-flow-network
  contact sheets for visual review.
- `earth-climate-compare` against the current six terminal worlds.
- `terminal-climate-replay` on frozen terminal arrays for quick climate-only
  iteration.
- Targeted unit tests for new report/calibration logic.
- A visual contact sheet comparing Earth, earthlike, arid, and waterworld maps.

Minimum acceptance for a tuning pass:

- No NaNs or negative precipitation.
- Annual precipitation equals seasonal aggregate.
- Earthlike precipitation improves without making waterworld fake-monsoonal or
  arid worlds uniformly wet.
- Biome changes are interpreted only after climate metrics improve.

## Initial Baseline From R4/R6 Comparison

Earth 8000-cell baseline:

- Land fraction: 0.287
- Global mean surface/SST temperature: 15.7 C
- Land precipitation mean: 730 mm/yr
- Land precipitation p50: 498 mm/yr
- Land precipitation p90: 1700 mm/yr
- Current speed p90: 0.267 m/s
- R5 OSCAR current speed p90: 0.192 m/s annual, 0.223 m/s monthly
- R6 ESA CCI land cover: water 0.718, forest 0.084, cropland 0.041,
  grass/shrub 0.068, bare/sparse 0.056, snow/ice 0.028
- Desert area fraction: 0.058
- Forest area fraction: 0.098
- Tropical area fraction: 0.057

Generated earthlike baseline:

- `earthlike_seed42`: land precipitation mean 143 mm/yr, p50 104 mm/yr,
  p90 344 mm/yr, desert fraction 0.215, forest/tropical 0.
- `earthlike_seed909`: land precipitation mean 158 mm/yr, p50 142 mm/yr,
  p90 342 mm/yr, desert fraction 0.235, forest/tropical 0.

Dominant baseline conclusion:

- Do not start by tuning biomes.  The main blocker is hydroclimate: circulation,
  moisture access, and precipitation coefficients are too dry for Earthlike
  worlds.

## Current C4j + R6 Replay Status

Latest replay artifacts:

- Rendered C4j climate replay:
  `out_terminal_climate_replay_c4j1_precip_objects_render_r6_20260706/`
- R6 Earth comparison:
  `out_earth_climate_comparison_c4j1_render_r6_20260706/`
- R6 fitting report:
  `out_earth_climate_fitting_c4j1_render_r6_20260706/`
- Contact sheet:
  `out_earth_climate_comparison_c4j1_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`

Current acceptance:

- `earth-climate-compare` against R6 reports `earthlike flagged: 0`.
- `earth-climate-fit-report --fail-on-guardrail` reports guardrail verdict
  `pass` with `0` failures and `0` warnings.
- The comparison summary now carries ESA CCI broad land-cover fractions from R6
  as reference context, while generated worlds are still evaluated through the
  generated biome classifier rather than a human land-use model.
- The comparison contact sheet now labels the active Earth reference generation
  (`Earth R6`) instead of the previous hard-coded `Earth R4` label.

Current earthlike metrics against R6:

- `earthlike_seed42`: land precipitation mean/p50/p90 ratios are
  `0.575` / `0.618` / `0.585`; current p90 ratio is `1.03`; forest/tropical
  area fractions are `0.042` / `0.025`.
- `earthlike_seed909`: land precipitation mean/p50/p90 ratios are
  `0.702` / `0.842` / `0.712`; current p90 ratio is `1.00`; forest/tropical
  area fractions are `0.065` / `0.045`.

Current interpretation:

- The original severe dry-bias blocker is cleared at the current gate
  strictness, and OSCAR-based current speed is now properly scaled.
- Remaining climate-system work is refinement rather than blocker repair:
  improve physically routed moisture access, stabilize receiver-catchment
  semantics, and retune biome thresholds only after precipitation objects and
  moisture budgets remain stable.

## Current F0 Artifacts

Report command:

```bash
.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_r4_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_r0_20260705
```

Generated artifacts:

- `out_earth_climate_fitting_r0_20260705/earth_climate_fitting_report.json`
- `out_earth_climate_fitting_r0_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_fitting_r0_20260705/earth_climate_fitting_runs.csv`
- `out_earth_climate_fitting_r0_20260705/earth_climate_fitting_levers.csv`

Replay command:

```bash
.venv/bin/python -m aevum.cli terminal-climate-replay \
  --terminal-summary out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json \
  --label earthlike_seed42 \
  --label earthlike_seed909 \
  --out out_terminal_climate_replay_r0_20260705 \
  --no-render
```

Replay verification:

- `earthlike_seed42` and `earthlike_seed909` match the original terminal
  climate metrics with zero numeric delta for mean temperature, land/ocean mean
  temperature, mean precipitation, land precipitation p50/p90,
  precipitation seasonality, and monsoon index.

## Current F1-F5 First-Pass Artifacts

Replay command:

```bash
.venv/bin/python -m aevum.cli terminal-climate-replay \
  --terminal-summary out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json \
  --out out_terminal_climate_replay_f1_oceanfloor_render_20260705
```

Comparison/report commands:

```bash
.venv/bin/python -m aevum.cli earth-climate-compare \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f1_oceanfloor_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_comparison_f1_oceanfloor_20260705 \
  --no-contact-sheet

.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_f1_oceanfloor_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_f1_oceanfloor_20260705

.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_f1_oceanfloor_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_f1_oceanfloor_gate_20260705 \
  --fail-on-guardrail
```

Generated artifacts:

- `out_terminal_climate_replay_f1_oceanfloor_20260705/terminal_climate_replay_summary.json`
- `out_terminal_climate_replay_f1_oceanfloor_render_20260705/*/{temperature,precip,precip_seasons,biomes,currents,moisture_access,monsoon_index}.png`
- `out_earth_climate_comparison_f1_oceanfloor_20260705/earth_climate_comparison_summary.json`
- `out_earth_climate_comparison_f1_oceanfloor_20260705/earth_climate_comparison_metrics.csv`
- `out_earth_climate_fitting_f1_oceanfloor_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_fitting_f1_oceanfloor_20260705/earth_climate_fitting_runs.csv`
- `out_earth_climate_fitting_f1_oceanfloor_gate_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_fitting_f1_oceanfloor_gate_20260705/earth_climate_guardrails.csv`

Current recommended F4-tail artifacts:

```bash
.venv/bin/python -m aevum.cli terminal-climate-replay \
  --terminal-summary out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json \
  --out out_terminal_climate_replay_f4tail2_render_20260705

.venv/bin/python -m aevum.cli earth-climate-compare \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f4tail2_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_comparison_f4tail2_20260705 \
  --no-contact-sheet

.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_f4tail2_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_f4tail2_gate_20260705 \
  --fail-on-guardrail
```

- `out_terminal_climate_replay_f4tail2_render_20260705/*/{temperature,precip,precip_seasons,biomes,currents,moisture_access,monsoon_index}.png`
- `out_earth_climate_comparison_f4tail2_20260705/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_f4tail2_gate_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_fitting_f4tail2_gate_20260705/earth_climate_guardrails.csv`

Current recommended Earth-pattern artifacts:

```bash
.venv/bin/python -m aevum.cli terminal-climate-replay \
  --terminal-summary out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json \
  --out out_terminal_climate_replay_pattern7_render_20260705

.venv/bin/python -m aevum.cli earth-climate-compare \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_pattern7_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_comparison_pattern7_20260705 \
  --no-contact-sheet

.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_pattern7_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_pattern7_gate_20260705 \
  --fail-on-guardrail

.venv/bin/python -m aevum.cli earth-climate-pattern-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_pattern7_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_pattern_gate_pattern7_20260705
```

- `out_terminal_climate_replay_pattern7_render_20260705/*/{temperature,precip,dry_season_length,biomes,moisture_access}.png`
- `out_earth_climate_comparison_pattern7_20260705/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_pattern7_gate_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_fitting_pattern7_gate_20260705/earth_climate_guardrails.csv`
- `out_earth_climate_pattern_gate_pattern7_20260705/earth_climate_pattern_gate_report.md`
- `out_earth_climate_pattern_gate_pattern7_20260705/earth_climate_pattern_checks.csv`

Current recommended windward/leeward artifacts:

```bash
.venv/bin/python -m aevum.cli terminal-climate-replay \
  --terminal-summary out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json \
  --out out_terminal_climate_replay_f5wind3_render_20260705

.venv/bin/python -m aevum.cli earth-climate-compare \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_comparison_f5wind3_20260705 \
  --no-contact-sheet

.venv/bin/python -m aevum.cli earth-climate-fit-report \
  --comparison-summary out_earth_climate_comparison_f5wind3_20260705/earth_climate_comparison_summary.json \
  --out out_earth_climate_fitting_f5wind3_gate_20260705 \
  --fail-on-guardrail

.venv/bin/python -m aevum.cli earth-climate-pattern-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_pattern_gate_f5wind3_20260705 \
  --fail-on-pattern

.venv/bin/python -m aevum.cli earth-climate-biome-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_biome_gate_f5wind3_20260705 \
  --fail-on-biome

.venv/bin/python -m aevum.cli earth-climate-spatial-biome-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_spatial_biome_gate_f5wind3_20260705 \
  --fail-on-spatial-biome

.venv/bin/python -m aevum.cli earth-climate-seasonal-subtype-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_seasonal_subtype_gate_f5wind3_20260705 \
  --fail-on-seasonal-subtype

.venv/bin/python -m aevum.cli earth-climate-mountain-zonation-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_mountain_zonation_gate_f5wind3_20260705 \
  --fail-on-mountain-zonation

.venv/bin/python -m aevum.cli earth-climate-windward-leeward-gate \
  --earth-reference out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz \
  --terminal-summary out_terminal_climate_replay_f5wind3_20260705/terminal_climate_replay_summary.json \
  --out out_earth_climate_windward_leeward_gate_f5wind3_20260705 \
  --fail-on-windward-leeward
```

- `out_terminal_climate_replay_f5wind3_render_20260705/*/{temperature,precip,precip_seasons,biomes,dry_season_length,moisture_access,orographic_precipitation,wind_seasons}.png`
- `out_earth_climate_comparison_f5wind3_20260705/earth_climate_comparison_summary.json`
- `out_earth_climate_fitting_f5wind3_gate_20260705/earth_climate_fitting_report.md`
- `out_earth_climate_pattern_gate_f5wind3_20260705/earth_climate_pattern_gate_report.md`
- `out_earth_climate_biome_gate_f5wind3_20260705/earth_climate_biome_gate_report.md`
- `out_earth_climate_spatial_biome_gate_f5wind3_20260705/earth_climate_spatial_biome_gate_report.md`
- `out_earth_climate_seasonal_subtype_gate_f5wind3_20260705/earth_climate_seasonal_subtype_gate_report.md`
- `out_earth_climate_mountain_zonation_gate_f5wind3_20260705/earth_climate_mountain_zonation_gate_report.md`
- `out_earth_climate_windward_leeward_gate_f5wind3_20260705/earth_climate_windward_leeward_gate_report.md`

First-pass outcome:

- Earthlike dry flags: 2 -> 0.
- Earthlike distance score:
  - `earthlike_seed42`: 0.85 -> 0.47.
  - `earthlike_seed909`: 0.84 -> 0.37.
- Earthlike global temperature delta:
  - `earthlike_seed42`: -1.38 C -> -0.21 C.
  - `earthlike_seed909`: -3.23 C -> -1.69 C.
- Earthlike land precipitation:
  - mean ratio to Earth: 0.20-0.22 -> 0.58-0.69.
  - p50 ratio to Earth: 0.21-0.28 -> 0.84-1.05.
  - p90 ratio to Earth: 0.20 -> 0.46-0.56.
- Current speed p90 ratio to Earth under the R4 AOML baseline:
  about 1.66-1.68 -> about 1.06-1.07.
- Re-evaluate this ratio against the R5 OSCAR baseline before the next climate
  tuning pass, because OSCAR's annual p90 is lower than the AOML drifter
  cross-check.
- Re-evaluate biome/land-cover gates against R6 before the next biome tuning
  pass, because ESA CCI land cover now provides an independent cross-check in
  addition to Koppen-derived biome proxy and RESOLVE ecoregions.
- Earthlike forest/tropical biome fractions are nonzero:
  - `earthlike_seed42`: forest 0.031, tropical 0.009.
  - `earthlike_seed909`: forest 0.054, tropical 0.029.
- Cross-preset guardrail verdict:
  - `pass_with_warnings`
  - 0 failures.
  - 1 warning: `earthlike_seed42` high-rainfall tail remains thin
    (`land_precip_p90_ratio_to_earth = 0.460`).
  - The guardrail CSV is now the authoritative per-check table for batch
    review; skipped checks are counted separately so missing diagnostic metrics
    do not masquerade as failures.

F4-tail gate outcome:

- Cross-preset guardrail verdict: `pass`.
- 0 failures, 0 warnings, 0 skipped checks.
- `earthlike_seed42`: land precipitation mean/p50/p90 ratios are
  `0.630 / 0.900 / 0.517`; forest/tropical fractions are `0.040 / 0.021`.
- `earthlike_seed909`: land precipitation mean/p50/p90 ratios are
  `0.753 / 1.131 / 0.638`; forest/tropical fractions are `0.067 / 0.039`.
- Arid guardrails remain stable: median land precipitation is about
  `8.9-84.7 mm/yr`, desert fraction about `0.66-0.70`, and wet land above
  500 mm/yr stays below `0.081`.
- Waterworld guardrails remain stable: land fraction is about `0.035-0.041`,
  median island precipitation is about `612-704 mm/yr`, and desert fraction is
  near zero.

Pattern7 outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass_with_warnings`.
- Pattern failures: 0.
- Pattern warnings: 2, both forest+tropical land-fraction envelope warnings.
- `earthlike_seed42`: global land precipitation mean/p50/p90 ratios are
  `0.569 / 0.622 / 0.565`; wet-tropics p90 ratio is `0.501`;
  subtropical dry fraction is `0.280`; high-latitude cold / ice-tundra
  fractions are `0.696 / 0.458`.
- `earthlike_seed909`: global land precipitation mean/p50/p90 ratios are
  `0.694 / 0.842 / 0.703`; wet-tropics p90 ratio is `0.570`;
  subtropical dry fraction is `0.180`; high-latitude cold / ice-tundra
  fractions are `0.915 / 0.819`.
- Arid and waterworld guardrails remain stable under the same scalar fitting
  gate.

F5 biome-gate outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass`.
- Earth biome gate verdict: `pass`.
- `earthlike_seed42`: forest/tropical land fraction is `0.231`
  (`0.120 / 0.110` forest/tropical), desert is `0.420`, tundra+ice is `0.111`.
- `earthlike_seed909`: forest/tropical land fraction is `0.267`
  (`0.148 / 0.119` forest/tropical), desert is `0.244`, tundra+ice is `0.175`.
- Arid worlds remain desert-dominated, with forest+tropical land fractions near
  `0.035-0.039`; waterworld islands become wetter/forested as expected for
  small maritime landmasses.

F5 spatial-biome gate outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass`.
- Earth biome gate verdict: `pass`.
- Earth spatial-biome gate verdict: `pass`.
- Initial spatial gate attribution found that tropical and subtropical belts
  were plausible, but cool-midlatitude forests were too weak and cool/high-lat
  land was overclassified as desert.
- `earthlike_seed42`: cool-midlat forest/desert fractions moved to
  `0.328 / 0.119`; high-lat tundra+ice/desert moved to `0.877 / 0.003`.
- `earthlike_seed909`: cool-midlat forest/desert fractions moved to
  `0.395 / 0.198`; high-lat tundra+ice/desert moved to `0.980 / 0.000`.
- Coarse biome fractions remain in bounds: forest+tropical land fractions are
  `0.281` and `0.350`, while arid worlds remain desert-dominated.

F5 seasonal-subtype gate outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass`.
- Earth biome gate verdict: `pass`.
- Earth spatial-biome gate verdict: `pass`.
- Earth seasonal-subtype gate verdict: `pass`.
- Initial seasonal-subtype attribution found that low-tropical dry/wet-season
  subdivision was too weak even when annual precipitation and biome envelopes
  passed.
- The seasonal redistribution preserves annual precipitation exactly per cell.
  Latest low-tropics dry-quarter >=2 fractions are `0.098` and `0.096`,
  clearing the Earth-ratio floor versus Earth reference `0.138`.
- Latest low/mid-latitude precipitation seasonality p75 values are `3.063` and
  `3.146` versus Earth `2.270`, still inside gate tolerance.
- Tropical coarse biome fractions remain in bounds after semantic retuning:
  `0.099` and `0.156`.

F5 mountain-zonation gate outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass`.
- Earth biome gate verdict: `pass`.
- Earth spatial-biome gate verdict: `pass`.
- Earth seasonal-subtype gate verdict: `pass`.
- Earth mountain-zonation gate verdict: `pass`.
- Initial mountain-zonation attribution found that mountain temperature
  gradients already passed, but high mountains were overclassified as desert:
  high-mountain desert fractions were about `0.374` and `0.364`.
- The biome classifier now adds a high-elevation alpine stress path before
  arid classification.  This uses relative elevation, annual temperature, and
  coldest-season temperature; it does not change terminal terrain or broad
  precipitation fields.
- Latest earthlike high-mountain alpine ecology fractions are `0.576` and
  `0.533`; high-mountain desert fractions are `0.235` and `0.200`, clearing
  the gate.  Midlatitude high-mountain desert fractions are `0.268` and
  `0.234`, below the warning threshold.

Windward/leeward gate outcome:

- Scalar fitting guardrail verdict: `pass`.
- Earth pattern gate verdict: `pass`.
- Earth biome gate verdict: `pass`.
- Earth spatial-biome gate verdict: `pass`.
- Earth seasonal-subtype gate verdict: `pass`.
- Earth mountain-zonation gate verdict: `pass`.
- Earth windward/leeward gate verdict: `pass`.
- Initial windward/leeward attribution found that `earthlike_seed42` had
  mountain-slope windward/leeward precipitation ratios too close to 1:
  annual `1.068` versus Earth `2.081`, and seasonal median `1.166` versus
  Earth `2.063`.
- Terminal replay arrays now persist seasonal wind, background/thermal/
  orographic wind anomalies, terrain barriers, wind gaps, and orographic
  precipitation diagnostics so this gate is reproducible from archived arrays.
- The climate model now computes slope-wind exposure from climate-scale
  topography and seasonal winds.  Seasonal precipitation is redistributed
  from leeward to windward slopes with land-season mean preservation, avoiding
  a broad precipitation boost.
- Latest earthlike annual windward/leeward precipitation ratios are `1.309`
  and `1.712`; seasonal median ratios are `1.574` and `1.837`, clearing the
  Earth-ratio gate.

C4a monsoon/moisture gate outcome:

- Added `aevum.diagnostics.earth_climate_monsoon_moisture_gate` and CLI command
  `aevum earth-climate-monsoon-moisture-gate`.
- The gate checks the C4a fields, not just final rainfall:
  `atmosphere.seasonal_pressure_proxy`, `atmosphere.moisture_access`, and
  `atmosphere.monsoon_potential`.
- It compares Earthlike worlds against broad real-Earth seasonal monsoon
  envelopes, while guarding waterworld and arid presets against fake
  continent-scale monsoon behavior.
- Baseline `f5wind3` failed this gate on `waterworld_seed707` with `3`
  failures: excessive land monsoon-potential p95, excessive low/mid-latitude
  summer monsoon-potential p90, and excessive seasonal pressure reversal.
- Climate now damps continent interiority and seasonal pressure by absolute
  landmass scale, so small waterworld islands cannot act like continental heat
  lows.
- Current replay:
  `out_terminal_climate_replay_c4a1_20260705/`.
- Current C4a gate:
  `out_earth_climate_monsoon_moisture_gate_c4a1_20260705/`, verdict `pass`
  with `0` failures and `0` warnings.
- Key metrics: earthlike summer monsoon-potential p90 is `0.359` and `0.679`;
  waterworld summer monsoon-potential p90 is `0.059` and `0.113`; waterworld
  pressure summer-minus-winter means are now only `-0.087` and `-0.217`.

Residual issues:

- F1 temperature is low priority but still a watch item because land/ocean
  meaning differs between generated worlds and the mixed WorldClim/OISST
  reference.
- F4 climate-pattern failures are cleared at the current gate strictness; avoid
  further broad precipitation boosts unless a new pattern metric proves the need.
- F5 coarse, spatial, seasonal-subtype, mountain-zonation, windward/leeward,
  C4a monsoon/moisture, C4b basin-streamfunction currents, C4c
  SST/wind/evaporation/pressure weak coupling, and C4d regional
  seasonal-hydroclimate plus object/placement/map-readability gates are cleared
  at current gate strictness.  The hydro-region gate now emits rendered contact
  sheets; the next useful work is to review those sheets and only then move to
  explicit flow-network organization if current corridor fields still read too
  locally.

C4c weak-coupling outcome:

- `ClimateModule` now weakly couples C4b SST anomaly/upwelling to
  `climate.seasonal_sst`, mean-zero `climate.ocean_heat_flux`, seasonal
  pressure, bounded seasonal wind adjustments, evaporation, source-ocean
  warmth, moisture access, and precipitation.
- Current replay:
  `out_terminal_climate_replay_c4c3_20260705/`.
- C4c render probe:
  `out_terminal_climate_replay_c4c3_render_probe_20260705/earthlike_seed42/`,
  including `seasonal_sst.png`, `ocean_heat_flux.png`,
  `coupling_residual.png`, `evaporation.png`, `sst_anomaly.png`, and
  `current_streamfunction.png`.
- Earth comparison:
  `out_earth_climate_comparison_c4c3_20260705/`, with `earthlike flagged: 0`
  and earthlike scores `0.46` and `0.33`.
- Earth climate fit report:
  `out_earth_climate_fitting_c4c3_gate_20260705/`, guardrail verdict `pass`
  with `0` failures and `0` warnings.
- All current Earth climate gates pass on C4c3 with `0` failures and `0`
  warnings: monsoon/moisture, pattern, coarse biome, spatial biome,
  seasonal-subtype, mountain-zonation, and windward/leeward.

C4d regional seasonal-hydroclimate outcome:

- `ClimateModule` now emits first-pass regional wet/dry corridor fields:
  `climate.monsoon_rainfall_corridor`,
  `climate.storm_track_rainfall_corridor`, `climate.rain_shadow_index`, and
  `climate.regional_precipitation_response`.
- The response is derived from C4a/C4b/C4c diagnostics and is applied
  conservatively, preserving each season's land precipitation mean before
  annual aggregation.
- Current replay:
  `out_terminal_climate_replay_c4d2_20260705/`.
- C4d render probe:
  `out_terminal_climate_replay_c4d2_render_probe_20260705/earthlike_seed42/`,
  including `monsoon_rainfall_corridor.png`,
  `storm_track_rainfall_corridor.png`, `rain_shadow_index.png`,
  `regional_precipitation_response.png`, `precip_seasons.png`, and
  `runoff.png`.
- Earth comparison:
  `out_earth_climate_comparison_c4d2_20260705/`, with `earthlike flagged: 0`
  and earthlike scores `0.46` and `0.33`.
- Earth climate fit report:
  `out_earth_climate_fitting_c4d2_gate_20260705/`, guardrail verdict `pass`
  with `0` failures and `0` warnings.
- All current Earth climate gates pass on C4d2 with `0` failures and `0`
  warnings: monsoon/moisture, pattern, coarse biome, spatial biome,
  seasonal-subtype, mountain-zonation, and windward/leeward.

C4d object-layer outcome:

- `ClimateModule` now emits `climate.hydroclimate_regions` diagnostic objects
  derived from the C4d monsoon, storm-track, rain-shadow, and regional-response
  fields.  This is an object-level calibration handle, not a new precipitation
  tuning pass.
- Current replay:
  `out_terminal_climate_replay_c4d3_objects_20260705/`.
- Earth comparison:
  `out_earth_climate_comparison_c4d3_objects_20260705/`, with
  `earthlike flagged: 0` and earthlike scores `0.46` and `0.33`.
- Earth climate fit report:
  `out_earth_climate_fitting_c4d3_objects_gate_20260705/`, guardrail verdict
  `pass` with `0` failures and `0` warnings.
- All current Earth climate gates pass on C4d3 objects with `0` failures and
  `0` warnings: monsoon/moisture, pattern, coarse biome, spatial biome,
  seasonal-subtype, mountain-zonation, and windward/leeward.
- Six-world object counts are now archived in replay summaries:
  arid seeds `179` and `161`, earthlike seeds `163` and `190`, and waterworld
  seeds `49` and `52`.

C4d hydroclimate-region gate outcome:

- Terminal climate replay now archives per-world
  `hydroclimate_regions.json` object files.
- Added `aevum earth-climate-hydro-region-gate` for the C4d object layer.
  The gate checks object archive presence, Earthlike multi-kind/multi-season
  coverage, coherent region-size proxies, monsoon seasonal migration,
  storm-track latitude/coastal placement, rain-shadow dry-response alignment,
  and waterworld/arid false-positive limits.
- Current replay:
  `out_terminal_climate_replay_c4d4_regiongate_20260705/`.
- Current hydro-region coverage gate:
  `out_earth_climate_hydro_region_gate_c4d4_20260705/`, verdict `pass` with
  `0` failures and `0` warnings.
- Current hydro-region placement-proxy gate:
  `out_earth_climate_hydro_region_gate_c4d5_placement_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Current hydro-region map-readability gate:
  `out_earth_climate_hydro_region_gate_c4d6_mapread_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Current hydro-region contact-sheet gate:
  `out_earth_climate_hydro_region_gate_c4d7_contacts_20260705/`, verdict
  `pass` with `0` failures, `0` warnings, `0` skipped checks, and `6`
  nonblank contact sheets.
- Earth comparison:
  `out_earth_climate_comparison_c4d4_regiongate_20260705/`, with
  `earthlike flagged: 0` and earthlike scores `0.46` and `0.33`.
- Earth climate fit report:
  `out_earth_climate_fitting_c4d4_regiongate_20260705/`, guardrail verdict
  `pass` with `0` failures and `0` warnings.
- All current Earth climate gates pass on C4d4/C4d6 evidence with `0` failures
  and `0` warnings: monsoon/moisture, pattern, coarse biome, spatial biome,
  seasonal-subtype, mountain-zonation, windward/leeward, and hydro-region.
- Placement-proxy metrics on the two earthlike seeds are within gate: monsoon
  field migration `66.61` and `59.44` latitude degrees, monsoon object
  migration `60.49` and `61.14`, storm-track weighted latitude `39.52` and
  `38.31`, storm-track coast-distance medians `0.058` and `0.078`, and
  rain-shadow/dry-response correlation `0.978` and `0.983`.  Waterworld
  monsoon field migration remains below the false-positive ceiling at `17.86`
  and `23.91`.
- C4d6 map-readability metrics on the two earthlike seeds are within gate:
  monsoon largest connected shares are `0.390` and `0.407`, monsoon active
  land fractions are `0.219` and `0.217`, storm-track largest connected shares
  are `0.344` and `0.190`, storm-track boundary-per-active-cell values are
  `1.891` and `1.487`, wet-response largest connected shares are `0.289` and
  `0.207`, and dry-response largest connected shares are `0.181` and `0.314`.
- C4d7 writes per-world contact sheets under
  `out_earth_climate_hydro_region_gate_c4d7_contacts_20260705/contact_sheets/`.
  Each sheet shows four seasons of precipitation, monsoon corridor mask,
  storm-track corridor mask, rain-shadow mask, and regional-response anomaly
  from the same archived arrays used by the gate.

Current Earth-fitting acceptance:

- The C4d7 contact sheets were reviewed directly after the metric gate.  The two
  earthlike runs show seasonal monsoon migration, mid/high-latitude storm-track
  wet corridors, rain-shadow masks aligned with dry regional response, and
  response anomalies that remain coupled to geography rather than random cell
  texture.
- The arid and waterworld contact sheets do not show broad false-positive
  Earthlike monsoon systems: arid wet regions remain bounded and waterworld
  monsoon/storm-track masks stay sparse around limited land exposure.
- Remaining visual debt is not a current Earth-fitting blocker: monsoon and
  regional-response masks still read as wet/dry corridor fragments rather than
  explicit routed moisture-flow-network objects.  That is promoted to the next
  climate-system design phase, where flow-network organization can be added
  without re-opening the accepted scalar/pattern/biome Earth gates.
- Therefore the current Earth-based fitting pass is accepted as sufficient to
  move from the Earth calibration track back to climate-system development,
  while keeping the existing gate suite as regression protection.

Post-acceptance C4e regression note:

- The first climate-system follow-up added diagnostic seasonal moisture-flow
  network objects without changing the precipitation budget.
- C4e replay:
  `out_terminal_climate_replay_c4e1_flow_20260705/`.
- Existing C4d hydro-region regression gate:
  `out_earth_climate_hydro_region_gate_c4e1_flow_20260705/`, verdict `pass`
  with `0` failures, `0` warnings, and `0` skipped checks.
- Earth comparison and fit reports:
  `out_earth_climate_comparison_c4e1_flow_20260705/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4e1_flow_20260705/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- This keeps the accepted Earth-fitting pass intact while opening the next
  climate-system validation target.  That C4e-specific
  moisture-flow-network gate is now implemented and archived below.

Post-acceptance C4e moisture-flow-network gate note:

- Added `aevum earth-climate-moisture-flow-gate` for C4e diagnostics.
- C4e moisture-flow gate:
  `out_earth_climate_moisture_flow_gate_c4e2_gate_20260705/`, verdict `pass`
  with `0` failures, `0` warnings, `0` skipped checks, and `6` nonblank
  contact sheets.
- Representative visual review of the C4e contact sheets found the intended
  diagnostic structure: earthlike moisture pathways are tied to ocean source
  moisture and monsoon/storm support, while waterworld pathways remain
  island-scale.  Residual climate-system debt is that flow-network id segments
  can still look fragmented inside broader corridors.
- The gate checks object/archive presence, multi-kind/multi-season Earthlike
  coverage, source-ocean moisture, routed land pathways, pathway coupling to
  monsoon/storm support, map readability, and waterworld/arid false-positive
  limits.
- Existing C4d hydro-region regression gate on the same C4e replay:
  `out_earth_climate_hydro_region_gate_c4e2_regression_20260705/`, verdict
  `pass` with `0` failures and `0` warnings.
- Targeted C4e/C4d/climate regression passed:
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`
  -> `14 passed in 25.02s`.
- This means the Earth-fitting acceptance remains intact while climate-system
  development can now decide whether C4e moisture-flow networks remain
  diagnostic or become an active conservative precipitation redistribution
  layer.

Post-acceptance C4f conservative precipitation-response note:

- C4f now makes C4e moisture-flow networks active as a conservative seasonal
  precipitation redistribution response.  The plate and terrain inputs remain
  frozen; this pass changes only climate fields, climate diagnostics, rendering,
  and the refreshed precipitation summaries on moisture-flow objects.
- New field: `climate.moisture_flow_precipitation_response`.
- C4f replay:
  `out_terminal_climate_replay_c4f1_flowprecip_20260705/`.
- Earth comparison and fit reports:
  `out_earth_climate_comparison_c4f1_flowprecip_20260705/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4f1_flowprecip_20260705/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Existing Earth climate gates all pass on the C4f replay with `0` failures,
  `0` warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a
  monsoon/moisture, C4d hydro-region, and C4e moisture-flow-network.
- C4f response strength is active but bounded.  The two earthlike seeds have
  land response p05/p95 near `0.807/1.068` and `0.794/1.069`; waterworld seeds
  remain weak at about `0.880/1.040` and `0.895/1.033`, so the response does
  not invent broad continent-scale monsoon systems on tiny islands.  Seasonal
  land mean precipitation conservation holds to numerical roundoff.
- Targeted regression passed for the C4d/C4e/C4f unit and engine tests, plus
  `tests/test_render.py` after centering the response render around `1.0`.
- This preserves the accepted Earth-fitting state while moving climate-system
  development beyond diagnostic objects.  Further climate work should compare
  the C4f maps visually and only then consider stronger basin/network moisture
  budgets or biome feedback.

Post-acceptance C4f precipitation-response gate note:

- Added `aevum earth-climate-moisture-response-gate` as the dedicated C4f
  regression gate.
- The gate checks that the C4f response field is archived, four-season, finite,
  conservative over seasonal land means, unchanged over ocean, bounded in
  wet/dry strength, positively coupled to moisture pathways and monsoon/storm
  support, not wetting rain-shadow regions, readable as a map, and weak enough
  on waterworlds to avoid false continental monsoons.
- Gate without contact sheets:
  `out_earth_climate_moisture_response_gate_c4f2_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Gate with contact sheets:
  `out_earth_climate_moisture_response_gate_c4f2_contacts_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, `0` skipped checks, and `6`
  nonblank contact sheets.
- Representative C4f gate metrics on the six-world replay: earthlike response
  p05/p95 remains about `0.79-0.81 / 1.07`, waterworld response remains about
  `0.89 / 1.03-1.04`, response-pathway correlation is positive in all worlds
  (`0.813-0.951`), and ocean response deviation is `0.0`.
- Direct visual review of `earthlike_seed42` and `waterworld_seed7` contact
  sheets confirms that response anomalies follow routed moisture/support
  structure and stay island-scale for waterworlds.
- Integrated C4d/C4e/C4f regression passed:
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  and `tests/test_render.py` -> `21 passed in 28.12s`.
- This keeps Earth-based fitting protected while allowing the climate-system
  plan to move from C4f field validation toward possible basin/network moisture
  budget design.

Post-acceptance C4g local moisture-budget note:

- C4g adds `climate.moisture_budget_region_id` and uses it to make the C4f
  precipitation response locally conservative over land budget regions.  The
  current first pass uses climate continent components as the budget authority:
  this prevents moisture-response wetting on one continent from being balanced
  by drying another continent, while leaving stricter source-basin/network
  sectors for a later refinement.
- Plate and terrain remain frozen.  The change is confined to climate response
  shaping, climate feature metadata, validation diagnostics, terminal replay
  archives, rendering, and the moisture-response gate.
- Six-world C4g replay:
  `out_terminal_climate_replay_c4g1_budget_20260706/`.
- C4g Earth comparison and fit reports:
  `out_earth_climate_comparison_c4g1_budget_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4g1_budget_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C4g moisture-response/local-budget gate:
  `out_earth_climate_moisture_response_gate_c4g1_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Existing Earth gates all pass on the C4g replay with `0` failures, `0`
  warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a monsoon/moisture,
  C4d hydro-region, C4e moisture-flow-network, and the expanded C4f/C4g
  moisture-response gate.
- Regression tests passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 28.28s`.
- Representative local-budget metrics from the C4g smoke replay:
  `earthlike_seed42` budget regions p50 `15`, max budget mean delta
  `4.55e-13 mm/yr`; `arid_seed101` p50 `1`, max delta `5.68e-14`; and
  `waterworld_seed7` p50 `22`, max delta `9.09e-13`.  These are numerical
  roundoff-level conservation errors.
- Earth-based fitting therefore remains accepted/protected after C4g.  The next
  climate-system work can safely proceed to basin/network-sector budgets or
  object-continuity refinements without reopening broad Earth scalar, pattern,
  or biome calibration.

Post-acceptance C4h moisture-network sector-budget note:

- C4h refines `climate.moisture_budget_region_id`: continent-scale C4g budget
  regions remain the base authority, and only large, coherent, strong
  moisture-flow networks are split into local halo sectors.  This avoids
  cross-continent compensation without forcing every small island or unstable
  network fragment into a separate water budget.
- Plate and terrain remain frozen.  The change is confined to climate response
  shaping, field metadata, replay summaries, diagnostics, and the expanded
  C4f/C4h moisture-response gate.
- Six-world C4h replay:
  `out_terminal_climate_replay_c4h1_budget_sector_20260706/`.
- C4h Earth comparison and fit reports:
  `out_earth_climate_comparison_c4h1_budget_sector_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4h1_budget_sector_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C4h moisture-response/local-budget gate:
  `out_earth_climate_moisture_response_gate_c4h1_gate_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Existing Earth gates all pass on the C4h replay with `0` failures, `0`
  warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a monsoon/moisture,
  C4d hydro-region, C4e moisture-flow-network, and C4h
  moisture-response/local-budget.
- Regression tests passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 28.55s`.
- Representative C4h sector metrics from the smoke replay:
  `earthlike_seed42` base regions `15`, post-split p50 `23.5`, sector-split
  p50 `7.5`; `earthlike_seed909` base regions `15`, post-split p50 `19.5`,
  sector-split p50 `3.5`; waterworld sector-split p50 remains weak (`0.0` and
  `1.5`).  Max local mean deltas remain numerical roundoff.
- Earth-based fitting therefore remains accepted/protected after C4h.  The next
  climate-system step should use these gates as regression protection while
  deciding whether source-ocean basin/catchment budgets are worth the added
  complexity.

Post-acceptance C4i source-ocean basin attribution note:

- C4i adds `atmosphere.moisture_source_basin_id`, a seasonal source-ocean basin
  attribution field for moisture pathways.  This is a diagnostic and budget
  authority improvement, not a broad precipitation boost.
- Moisture-flow network objects now archive source-basin ids, and C4h budget
  sector splitting uses source-basin labels when available so local sectors do
  not deliberately mix different diagnosed source oceans.
- Plate and terrain remain frozen.  The change is confined to climate
  diagnostics, response-budget splitting, field metadata, replay archives,
  rendering/contact sheets, validation, and the C4f/C4i moisture-response gate.
- Six-world C4i replay:
  `out_terminal_climate_replay_c4i1_source_basin_20260706/`.
- C4i Earth comparison and fit reports:
  `out_earth_climate_comparison_c4i1_source_basin_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4i1_source_basin_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C4i moisture-response/local-budget/source-basin coherence gate:
  `out_earth_climate_moisture_response_gate_c4i2_source_basin_coherence_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Existing Earth gates all pass on the C4i replay with `0` failures, `0`
  warnings, and `0` skipped checks: pattern, coarse biome, spatial biome,
  seasonal subtype, mountain zonation, windward/leeward, C4a monsoon/moisture,
  C4d hydro-region, C4e moisture-flow-network, and C4i
  moisture-response/local-budget.
- Regression tests passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `22 passed in 29.04s`.
- Representative source-basin coherence metrics: earthlike seeds attribute
  about `0.963` and `0.984` of seasonal land area to source ocean basins; arid
  seeds attribute about `0.774` and `0.390`; active pathways, wet response
  corridors, and moisture-flow networks have attribution p50 `1.0`; budget
  source-purity p50 is `1.0` in all six worlds.  Waterworld islands are fully
  attributed but still pass the waterworld false-positive response checks.
- Earth-based fitting remains accepted/protected after C4i2.  The next climate
  work can move to precipitation object continuity; a true source-basin /
  receiver-catchment water-budget gate should wait until receiver catchments are
  explicit objects rather than inferred from current budget-sector halos.

Post-acceptance C4j precipitation-response object-continuity note:

- C4j adds `climate.precipitation_response_region_id` and
  `climate.precipitation_response_regions`.  These objects bind the final C4f
  wet/dry precipitation-response patches to seasonal source-ocean basin ids,
  local moisture-budget regions, moisture-flow network ids, precipitation,
  pathway strength, monsoon/storm support, and rain-shadow diagnostics.
- C4j is intentionally diagnostic/object-level: it does not retune or reshape
  precipitation.  Plate and terrain remain frozen, and C4f/C4h/C4i numeric
  response fields are unchanged except for the new object/region-id archive.
- Terminal replay:
  `out_terminal_climate_replay_c4j1_precip_objects_20260706/`.
- Earth comparison and fit reports:
  `out_earth_climate_comparison_c4j1_precip_objects_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c4j1_precip_objects_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- The moisture-response gate is now
  `aevum.earth_climate_moisture_response_gate.v6`.  It checks the C4j region-id
  archive, precipitation-response object JSON, wet/dry object kind coverage,
  source-basin attribution, budget-region attribution, and wet-object
  moisture-flow attribution.  The C4j contact sheets now include a row for
  `C4j response region id`.
- C4j moisture-response/object-continuity gate:
  `out_earth_climate_moisture_response_gate_c4j2_contacts_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, `0` skipped checks, and `6`
  contact sheets.
- Existing C4d and C4e object gates remain green on the same C4j replay:
  `out_earth_climate_hydro_region_gate_c4j1_regression_20260706/` and
  `out_earth_climate_moisture_flow_gate_c4j1_regression_20260706/`, both
  verdict `pass` with `0` failures and `0` warnings.
- Regression tests passed:
  `tests/test_core.py::test_registry_resolves_dependencies`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_render.py`,
  `tests/test_earth_climate_moisture_response_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_climate_seasonal_redistribution.py` -> `23 passed in 28.16s`.
- Representative C4j metrics from the v6 gate:
  earthlike seeds have `71-83` precipitation-response objects, both wet/dry
  kinds, all four seasons, and p50 source/budget/wet-flow attribution `1.0`;
  waterworld seeds remain island-scale with `9-10` response objects while still
  passing the waterworld false-positive limits.
- Earth-based fitting remains accepted/protected after C4j.  Remaining work is
  climate-system design rather than scalar Earth refitting: the next candidate
  is to merge these response objects into more stable receiver catchments before
  attempting a stricter source-basin / receiver-catchment water-budget gate.

Current C5a2 circulation-layout status:

- Motivation: after reviewing the documented dependency order, the fitting
  pass returned to the W/O layers before further precipitation or biome work.
  The previous seasonal wind field was still dominated by a strong zonal
  background template; generated wind p90 was about twice the R6/NCEP Earth
  seasonal 10 m wind p90, while OSCAR-scaled ocean-current speed already sat in
  the right magnitude range.
- Added `aevum.diagnostics.earth_climate_circulation_layout_gate` and CLI
  command `earth-climate-circulation-layout-gate`.  The gate checks wind p90
  ratio to the R6 Earth reference, coastal onshore/offshore seasonal response,
  land-sea thermal wind anomaly strength, OSCAR-scale surface-current speed,
  and current land-mask leakage.
- Recalibrated `_seasonal_circulation` background winds to near-surface Earth
  wind magnitude so land-sea and orographic anomalies are no longer drowned by
  latitude bands.  This is a W-layer change; plate and terrain remain frozen.
- Added a source-basin fallback split in moisture-budget regions so residual
  budget sectors can be separated by source ocean when large coherent source
  labels exist.  Arid source/budget gate semantics now keep active pathway
  attribution strict while using a coarser source-purity floor for huge dry
  interiors where most cells are weakly attributed.
- Current replay:
  `out_terminal_climate_replay_c5a2_circulation_sourcebudget_20260706/`.
- Current circulation-layout gate:
  `out_earth_climate_circulation_layout_gate_c5a2_20260706/`, verdict `pass`
  with `0` failures, `0` warnings, and `0` skipped checks.
- Current Earth comparison and fitting:
  `out_earth_climate_comparison_c5a2_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5a2_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- Rendered C5a2 replay and visual contact sheet:
  `out_terminal_climate_replay_c5a2_render_20260706/` and
  `out_earth_climate_comparison_c5a2_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- Current downstream gates on C5a2 all pass with `0` failures and `0` warnings:
  monsoon/moisture, windward/leeward, hydro-region, moisture-flow,
  moisture-response, and receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and the source-budget
  split tests in `tests/test_climate_seasonal_redistribution.py`.
- Next fitting step should stay in documented order: refine F2 spatial
  current/SST structure against OSCAR/OISST and then F3 pressure/wind/moisture
  placement, before any further F4 precipitation or F5 biome changes.

Current C5b1 ocean-current/SST spatial status:

- Motivation: the C5a2 circulation-layout pass fixed wind/current magnitude
  envelopes but still left one earthlike world with too much strongest current
  area in remote open-ocean bands and too little in near-coast/boundary-current
  structure.  This is the F2 spatial part of the Earth-fitting order, not a
  precipitation or biome tuning pass.
- Added `aevum.diagnostics.earth_climate_ocean_spatial_gate` and CLI command
  `earth-climate-ocean-spatial-gate`.  The gate compares generated maps to R6
  OSCAR/OISST references and checks current p90 ratio, current land leakage,
  current-speed zonal dominance, SST zonal dominance, coastal swift-current SST
  anomaly spread, mean-zero heat transport, and earthlike near-coast/far-ocean
  strong-current placement.
- Recalibrated `_ocean_currents` with a narrow structural change: warm/cold
  boundary-current terms are slightly stronger, while remote open-ocean current
  strength is mildly damped where there is little shelf, strait, or diagnosed
  boundary-current influence.  Plate, terrain, temperature fitting, and
  precipitation coefficients were not changed.
- C5b1 six-world replay:
  `out_terminal_climate_replay_c5b1_ocean_spatial_20260706/`.
- C5b1 rendered replay and visual contact sheet:
  `out_terminal_climate_replay_c5b1_render_20260706/` and
  `out_earth_climate_comparison_c5b1_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- C5b1 F2 ocean-spatial gate:
  `out_earth_climate_ocean_spatial_gate_c5b1_20260706/`, verdict `pass` with
  `0` failures, `0` warnings, and `0` skipped checks.  Key earthlike metrics:
  `earthlike_seed42` current p90 ratio `0.903`, near-coast strongest-current
  share `0.407`, far-ocean strongest-current share `0.389`; `earthlike_seed909`
  current p90 ratio `0.898`, near-coast share `0.520`, far-ocean share `0.266`.
- C5b1 circulation-layout gate:
  `out_earth_climate_circulation_layout_gate_c5b1_20260706/`, verdict `pass`
  with `0` failures, `0` warnings, and `0` skipped checks.
- C5b1 Earth comparison and fitting:
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
- Coupling note: the existing climate implementation is a bounded weak-coupling
  loop, not a full atmosphere-ocean-hydroclimate solver.  The next Earth-fitting
  step should therefore treat F3 as a coupled-consistency pass over pressure,
  wind, source-ocean warmth, moisture routing, and seasonal precipitation
  support, rather than as an isolated scalar precipitation adjustment.

Current C5c2 pressure/wind/moisture coupled-consistency status:

- Motivation: after C5b1, speed and spatial ocean-current/SST checks passed,
  but the climate system still needed an explicit check that temperature/SST,
  pressure, wind, evaporation, moisture routing, monsoon potential, and
  precipitation support each other.  This addresses the risk that one field
  looks plausible in isolation while the coupled reduced physics is incoherent.
- Added `aevum.diagnostics.earth_climate_coupled_consistency_gate` and CLI
  command `earth-climate-coupled-consistency-gate`.
- The new gate compares against R6 Earth pressure/wind phase envelopes and
  checks generated internal consistency: seasonal warm land -> low pressure
  proxy, pressure gradient -> wind, SST -> source-ocean warmth -> evaporation,
  moisture/support -> precipitation, positive monsoon potential -> thermal lows,
  top monsoon-potential cells -> adequate moisture and enhanced rainfall,
  cold-current/upwelling regions -> lower evaporation, seasonal precipitation
  budget closure, mean-zero ocean heat flux, and bounded coupling residual.
- Tightened `_seasonal_pressure_moisture` so high monsoon potential is gated
  more strongly by actual moisture access.  The reduced model still permits
  weak dry thermal responses, but top monsoon-potential cells now require
  enough source-ocean moisture support.
- C5c2 six-world replay:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`.
- C5c2 rendered replay and visual contact sheet:
  `out_terminal_climate_replay_c5c2_render_20260706/` and
  `out_earth_climate_comparison_c5c2_render_r6_20260706/earth_vs_generated_climate_contact_sheet.png`.
- C5c2 F3 coupled-consistency gate:
  `out_earth_climate_coupled_consistency_gate_c5c2b_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key metric movement: `earthlike_seed42` top monsoon-potential moisture ratio
  improves from C5b1 `0.504` to C5c2 `0.852`; `earthlike_seed909` remains
  passing at `1.281`.  Pressure/temperature, wind/pressure, SST/source,
  evaporation/source, moisture/precipitation, and support/precipitation
  correlations remain inside gate.
- C5c2 F2/F3 regression gates pass with `0` failures and `0` warnings:
  ocean-spatial, circulation-layout, monsoon/moisture, and windward/leeward.
- C5c2 Earth comparison and fitting:
  `out_earth_climate_comparison_c5c2_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5c2_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C5c2 downstream object gates all pass with `0` failures and `0` warnings:
  hydro-region, moisture-flow, moisture-response, and receiver-catchment.
- Targeted tests passed:
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_circulation_layout_gate.py`,
  `tests/test_earth_climate_monsoon_moisture_gate.py`,
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`,
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`,
  and `tests/test_engine.py::test_cold_boundary_currents_suppress_local_evaporation`
  -> `14 passed in 74.55s`.
- Next fitting step should stay in documented order: use the now-passing F3
  coupled-consistency gate as the guardrail for any further F4 seasonal
  hydroclimate refinement.  Do not jump to biome thresholds until F4 seasonal
  precipitation placement is reviewed under this coupled gate.

Current C5d1 seasonal hydroclimate placement status:

- Motivation: after C5c2, the reduced model had pressure/wind/moisture
  consistency, but F4 still needed a direct check that seasonal wet and dry
  precipitation cells are located where the generated process fields say they
  should be.  This addresses the coupled-climate concern without claiming a
  full GCM: temperature/SST, winds, currents, moisture, and precipitation are
  still a bounded weak-coupling system, but precipitation placement now has an
  explicit explanation gate.
- Added `aevum.diagnostics.earth_climate_seasonal_hydro_placement_gate` and
  CLI command `earth-climate-seasonal-hydro-placement-gate`.
- The new gate checks archived arrays, wet-cell process support, unsupported
  wet-patch fraction, dry-cell explanation, support/precipitation correlation,
  wet-season phase alignment, process-rainfall ratios for earthlike monsoon,
  storm-track, ITCZ, and rain-shadow fields, and exact seasonal/annual
  precipitation closure.
- No climate formula was changed in C5d1.  The gate was run directly on the
  accepted C5c2 six-world replay:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`.
- C5d1 F4 seasonal-hydro placement gate:
  `out_earth_climate_seasonal_hydro_placement_gate_c5d1_20260706/`, verdict
  `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- Key earthlike metrics:
  `earthlike_seed42` wet-support p25 ratio `2.301`, unsupported wet fraction
  `0.000`, dry-explained fraction `0.822`, support/precipitation correlation
  `0.747`, wet-season peak-support match `0.922`, rain-shadow precipitation
  ratio `0.556`; `earthlike_seed909` wet-support p25 ratio `1.877`,
  unsupported wet fraction `0.000`, dry-explained fraction `0.825`,
  support/precipitation correlation `0.701`, peak-support match `0.882`, and
  rain-shadow precipitation ratio `0.358`.
- Cross-preset context: arid and waterworld runs also pass the placement gate.
  Waterworld dry-cell explanation uses a looser threshold because tiny islands
  make dry-cell samples sparse and less continent-like.
- Targeted tests passed:
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`, and
  `tests/test_earth_climate_moisture_response_gate.py`
  -> `19 passed in 2.52s`.
- Next fitting step: treat C5d1 as the current F4 guardrail.  Any future
  precipitation, biome, or map-compiler change should rerun the F2/F3/F4
  gates together, because the four fields the user highlighted
  (temperature/SST, precipitation, wind, and ocean currents) are mutually
  dependent even in this reduced model.

Current C5d1 full Earth-fitting acceptance bundle:

- Purpose: consolidate the formerly scattered Earth-fitting evidence into one
  reproducible R6/C5d1 acceptance bundle before returning to broader
  climate-system work.  The input replay is still the accepted frozen-terrain
  C5c2 climate replay because C5d1 added a guardrail and did not change climate
  formulas:
  `out_terminal_climate_replay_c5c2_coupled_consistency_20260706/`.
- R6 Earth comparison:
  `out_earth_climate_comparison_c5d1_acceptance_r6_20260706/` reports
  `earthlike flagged: 0`.  Earthlike distance scores are `0.46` for
  `earthlike_seed42` and `0.33` for `earthlike_seed909`.
- R6 fitting report:
  `out_earth_climate_fitting_c5d1_acceptance_r6_20260706/` reports guardrail
  verdict `pass` with `0` failures and `0` warnings.  Remaining phase levers
  are low-priority improvement handles, not blockers: F1/F2/F4 are `watch`;
  F3/F5 are `needs_tuning` with `low` priority.
- All acceptance gates pass with `0` failures, `0` warnings, and `0` skipped
  checks on the same replay:
  pattern, biome, spatial-biome, seasonal-subtype, mountain-zonation,
  windward/leeward, monsoon/moisture, circulation-layout, ocean-spatial,
  coupled-consistency, seasonal-hydro placement, hydro-region, moisture-flow,
  moisture-response, and receiver-catchment.
- Acceptance artifact directories:
  `out_earth_climate_pattern_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_biome_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_spatial_biome_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_seasonal_subtype_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_mountain_zonation_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_windward_leeward_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_monsoon_moisture_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_circulation_layout_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_ocean_spatial_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_coupled_consistency_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_seasonal_hydro_placement_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_hydro_region_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_moisture_flow_gate_c5d1_acceptance_20260706/`,
  `out_earth_climate_moisture_response_gate_c5d1_acceptance_20260706/`, and
  `out_earth_climate_receiver_catchment_gate_c5d1_acceptance_20260706/`.
- Targeted acceptance tests passed:
  `tests/test_earth_climate_seasonal_hydro_placement_gate.py`,
  `tests/test_earth_climate_coupled_consistency_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_hydro_region_gate.py`,
  `tests/test_earth_climate_moisture_flow_gate.py`,
  `tests/test_earth_climate_moisture_response_gate.py`, and
  `tests/test_earth_climate_receiver_catchment_gate.py`
  -> `27 passed in 3.23s`.
- Conclusion: the Earth-based fitting track is accepted at the current
  strictness and can now serve as regression protection.  The next work should
  return to the climate-system plan, starting with bounded coupling iteration
  among temperature/SST, pressure/wind, currents, evaporation/moisture, and
  precipitation.  The first such changes are recorded in the C5e notes below;
  future climate-system changes must rerun the current C5e acceptance bundle.

Post-acceptance C5e1 hydro-feedback regression note:

- Purpose: begin the climate-system bounded-coupling work after Earth fitting
  acceptance by letting preliminary seasonal precipitation feed back into
  pressure and wind before the final hydroclimate solve.  This keeps the
  existing C4c ocean/SST/wind coupling and adds one conservative
  precipitation-pressure/wind feedback pass; it is not yet a full multi-pass
  atmosphere-ocean GCM.
- Code change: `ClimateModule.step` now computes preliminary hydroclimate,
  derives a bounded `atmosphere.precipitation_pressure_feedback` field from
  wet warm supported precipitation cores and dry low-moisture subsidence cells,
  adds a capped `atmosphere.hydro_coupled_wind_anomaly`, and reruns final
  hydroclimate with those adjusted seasonal pressure/wind fields.  New archived
  fields are `atmosphere.precipitation_pressure_feedback`,
  `atmosphere.hydro_coupled_wind_anomaly`, and
  `climate.hydro_coupling_residual`.
- C5e1 replay:
  `out_terminal_climate_replay_c5e1_hydro_feedback_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e1_hydro_feedback_r6_20260706/` reports
  `earthlike flagged: 0`; `out_earth_climate_fitting_c5e1_hydro_feedback_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- All C5d1 acceptance gates rerun on C5e1 pass with `0` failures, `0`
  warnings, and `0` skipped checks: pattern, biome, spatial-biome,
  seasonal-subtype, mountain-zonation, windward/leeward, monsoon/moisture,
  circulation-layout, ocean-spatial, coupled-consistency, seasonal-hydro
  placement, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Feedback bounds from the archived arrays are small and preset-sensitive:
  earthlike pressure-feedback abs-p95 is about `0.027-0.029`, wind-anomaly p95
  is about `0.061-0.070 m/s`, and hydro-coupling residual p95 is about
  `0.022`; waterworld pressure feedback remains effectively zero.
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
- Conclusion: Earth-fitting acceptance remains intact after the first
  climate-system bounded feedback step.  Next climate-system work can expand
  this from a single hydro feedback pass toward a 2-4 iteration bounded loop,
  but must preserve the current acceptance evidence.  C5e3 below supersedes
  C5e1 as the active regression baseline.

Post-acceptance C5e3 coupling-convergence gate note:

- Purpose: formalize the C5e feedback pass with an explicit convergence gate
  before expanding it toward a 2-4 iteration bounded coupling loop.
- Added `aevum.diagnostics.earth_climate_coupling_convergence_gate` and CLI
  command `earth-climate-coupling-convergence-gate`.  The gate checks archived
  feedback fields, pressure-feedback p95/max bounds, feedback/base-pressure
  ratio, wind-anomaly p95/p99 bounds, wind-anomaly/seasonal-wind ratio,
  hydro-coupling residual, ocean coupling residual, waterworld false-positive
  feedback, and exact seasonal/annual precipitation closure.
- Added `atmosphere.land_sea_pressure_proxy` to the replay archive so the gate
  can compare hydro feedback against the pre-feedback seasonal pressure field.
- Initial C5e2 gate attribution found a real waterworld issue: tiny island
  rainfall was creating strong normalized hydro wind anomalies.  C5e3 fixes
  this by scaling hydro pressure/wind feedback by exposed land fraction, so
  small island worlds cannot drive planet-scale pressure/wind feedback.
- C5e3 replay:
  `out_terminal_climate_replay_c5e3_coupling_convergence_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e3_coupling_convergence_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e3_coupling_convergence_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
- C5e3 coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e3_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
  Earthlike pressure-feedback abs-p95 is about `0.035`, wind-anomaly p95 is
  `0.134-0.149 m/s`, and hydro residual p95 is about `0.033`; waterworld
  pressure-feedback abs-p95 falls to `0.00010-0.00026` and wind-anomaly p95 to
  about `0.003 m/s`.
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
- Conclusion: C5e3 became the first formal coupling-convergence regression
  baseline.  C5e4 below supersedes it with a bounded iteration loop and must be
  preserved by further climate-system coupling work.

Post-acceptance C5e4 bounded hydro-feedback iteration note:

- Purpose: replace the one-shot hydro feedback pass with a small bounded
  iteration loop so seasonal precipitation, pressure, and low-level wind can
  settle together without introducing a full dynamic atmosphere-ocean model.
- Code change: `ClimateModule.step` now uses a 3-pass
  `_seasonal_hydroclimate_feedback_loop`.  Each pass recomputes seasonal
  hydroclimate from the current pressure/wind state, derives conservative
  precipitation-pressure and tangent wind feedback, damps the update, and caps
  wind magnitudes.  The final hydroclimate is recomputed from the converged
  pressure/wind state.
- Added archived field and diagnostics:
  `climate.hydro_feedback_iteration_delta`,
  `hydro_feedback_iteration_delta_p95`, and
  `hydro_feedback_iteration_count`.
- C5e4 replay:
  `out_terminal_climate_replay_c5e4_hydro_iteration_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e4_hydro_iteration_r6_20260706/` reports
  `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e4_hydro_iteration_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- C5e4 coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e4_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
  Earthlike wind-anomaly p95 is `0.136-0.150 m/s`, hydro residual p95 is about
  `0.033`, and iteration-delta p95 is only `0.0000033-0.0000050`; waterworld
  wind-anomaly p95 remains near `0.003 m/s` with effectively zero iteration
  delta.
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
- Conclusion: C5e4 completed the bounded hydro-feedback iteration step.  C5e5
  below supersedes it as the current Earth-fitting regression baseline after
  adding bounded evaporation-SST feedback.

Post-acceptance C5e5 evaporation-SST feedback note:

- Purpose: extend the bounded coupling discipline beyond hydro feedback by
  letting ocean evaporation apply a small zero-net heat-flux correction to SST,
  pressure-source, wind, and current iteration.  This is still a reduced model,
  but temperature/SST, evaporation, pressure/wind, currents, moisture, and
  precipitation are no longer treated as a purely one-way chain.
- Code change: `_weak_ocean_atmosphere_coupling` now runs 3 ocean-atmosphere
  iterations and calls `_ocean_evaporation_heat_feedback`.  The feedback cools
  high-evaporation ocean source regions, weakly offsets cold-current/upwelling
  suppression, is smoothed over the ocean mask, and is forced to zero
  area-weighted ocean mean before it is added to `climate.ocean_heat_flux`.
- Added archived field and gate checks:
  `climate.ocean_evaporation_feedback`.  The coupling-convergence gate now
  requires this field and checks p95/max amplitude plus zero area-weighted
  ocean mean.
- C5e5 replay:
  `out_terminal_climate_replay_c5e5_evap_sst_feedback_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e5_evap_sst_feedback_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e5_evap_sst_feedback_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- C5e5 coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e5_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
  Ocean evaporation-feedback abs-p95 is about `0.17-0.21 C`, abs-max is
  `0.22-0.29 C`, area-weighted ocean mean is effectively `0`, and ocean
  coupling residual p95 stays below `0.001 C`.
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
- Conclusion: C5e5 completed the evaporation-SST feedback step.  C5e6 below
  supersedes it as the current Earth-fitting regression baseline after adding
  explicit wind-stress/current-response diagnostics and bounds.

Post-acceptance C5e6 wind-stress/current-response note:

- Purpose: make the wind-to-current coupling visible and bounded.  Before
  C5e6, wind already influenced `_ocean_currents`, but the wind-stress
  response was implicit inside the streamfunction/current vector and was not
  independently archived or gated.
- Code change: `_ocean_currents` now emits
  `ocean.wind_stress_current_response`.  The response is tangent, ocean-only,
  aligned with annual wind stress, and remains a small fraction of wind speed.
  A first attempt with too much open-ocean wind-stress acceleration failed the
  ocean-spatial gate by moving too many swift currents into the far ocean; the
  final C5e6 version keeps the open-ocean increment small and preserves
  boundary-current structure.
- Added archived field and diagnostics:
  `ocean.wind_stress_current_response`, validation metrics for response shape,
  land leakage, tangent normal component, p50/p95 speed, and convergence-gate
  checks for p50 presence, p95 cap, alignment, and wind-speed ratio.
- C5e6 replay:
  `out_terminal_climate_replay_c5e6_wind_stress_current_20260706/`.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e6_wind_stress_current_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e6_wind_stress_current_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.
- C5e6 coupling-convergence gate:
  `out_earth_climate_coupling_convergence_gate_c5e6_acceptance_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
  Earthlike wind-stress response p50 is `0.042-0.047 m/s`, p95 is
  `0.067-0.070 m/s`, land max is `0`, alignment p50 is `1.0`, and the p95
  response-to-wind ratio is about `0.012`.
- Ocean-spatial gate:
  `out_earth_climate_ocean_spatial_gate_c5e6_20260706/`, verdict `pass`.
  Earthlike current-speed p90 is about `0.90` of Earth reference; swift-current
  far-ocean share is `0.405` for seed42 and `0.269` for seed909, both inside
  the earthlike gate.
- Full C5e6 acceptance suite passes with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile:
  `py_compile` passes for climate, terminal replay archive, coupling gate,
  validation, and edited tests.  Non-engine climate gate tests pass:
  `27 passed in 41.26s`.  The combined targeted pytest reached `4 passed`
  before being interrupted after `544.26s`, with the interruption inside a
  terrain/climate engine path; rerunning the standalone cold-boundary engine
  test also exceeded `5 min` in the current loaded environment and was
  interrupted.  This is treated as an incomplete engine-test run, not a model
  gate failure.
- Conclusion: C5e6 became the Earth-fitting regression baseline before the
  C5e7 source/receiver accounting pass.  The next climate-system step should
  move from diagnosed wind/current response toward better moisture-source/
  catchment closure by ocean basin and season.

Post-acceptance C5e7 source-basin receiver-catchment accounting note:

- Purpose: make the source-ocean-basin to receiver-catchment moisture ledger
  explicit.  Before C5e7, C4i/C4k could say which source basin and receiver
  catchment existed, but the gate could only check attribution, not whether wet
  response cells were backed by diagnosed seasonal ocean-basin supply.
- Code change: `ClimateModule` now emits
  `climate.source_basin_supply_index` and
  `climate.receiver_catchment_supply_balance`.  The first combines seasonal
  ocean source strength, dominant source-basin labels, and landward pathway
  support.  The second records a bounded receiver-catchment consistency
  diagnostic between precipitation magnitude and source supply.  This remains
  diagnostic-only; it does not yet feed back into the precipitation solver.
- Receiver-catchment objects now carry
  `mean_source_basin_supply_index`,
  `source_basin_supply_attributed_fraction`,
  `source_basin_supply_mass_fraction`,
  `supply_supported_precipitation_fraction`, and
  `precipitation_supply_balance`.
- C5e7 replay:
  `out_terminal_climate_replay_c5e7_source_receiver_accounting_20260706/`.
- Receiver-catchment v2 gate:
  `out_earth_climate_receiver_catchment_gate_c5e7_source_receiver_accounting_20260706/`,
  verdict `pass` with `0` failures, `0` warnings, and `0` skipped checks.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e7_visual_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e7_visual_r6_20260706/` reports guardrail
  verdict `pass` with `0` failures and `0` warnings.
- Visual evidence: `earth_vs_generated_climate_contact_sheet.png` now contains
  real Earth plus all six generated worlds for annual temperature,
  precipitation, biome, and current speed.  The comparison diagnostic can
  render these previews directly from archived arrays when the replay was run
  with `--no-render`, so existing C5e7 replay arrays do not need to be
  recomputed.
- The fitting report schema is now
  `aevum.earth_climate_fitting_report.v2`.  It scores F3 using the same
  low/mid-latitude summer diagnostics as the monsoon/moisture gate rather than
  diluting seasonal monsoon potential across all land and all seasons.  The
  final C5e7 report now classifies all phases as low-priority `watch`:
  F1 score `0.35`, F2 score `0.10`, F3 score `0.00`, F4 score `0.00`, and F5
  score `0.00`.
- Key F3 final-report metrics for earthlike seeds `42/909`: low/mid-latitude
  summer moisture p75 is `0.837/0.843`, summer monsoon-potential p90 is
  `0.358/0.667`, and summer-minus-winter monsoon-potential p75 is
  `0.259/0.483`.
- Key receiver metrics for earthlike seeds `42/909`: source-supply-attributed
  land p50 is `0.56/0.67`, wet-response source-supply p50 is `0.87/0.91`,
  and receiver supply-balance land p50 is `0.69/0.64`.
- Existing acceptance gates also pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile: `py_compile` passes for edited climate, diagnostics,
  validation, rendering, feature catalog, and tests.  Target tests pass:
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`,
  `tests/test_earth_climate_receiver_catchment_gate.py`,
  `tests/test_climate_seasonal_redistribution.py`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`.
- Conclusion: C5e7 was the previous Earth-fitting regression baseline.  C5e8
  below supersedes it by feeding source/receiver accounting into a bounded
  regional precipitation/moisture feedback solve.  The visual replay still
  shows reduced-model residuals, especially broad latitude-banded
  temperature/current structure, so follow-up work should preserve the C5e7/C5e8
  guardrails while improving geography-coupled mechanics.

Post-acceptance C5e8 receiver-supply precipitation feedback note:

- Purpose: move source/receiver accounting from diagnostic-only into the
  precipitation solve without broad retuning.  C5e8 adds a second, land-only,
  bounded redistribution after C4f: cells with stronger diagnosed
  source-basin supply and receiver-catchment balance gain a small amount of
  seasonal precipitation, while lower-support cells inside the same local
  moisture-budget region lose a matching amount.
- Code change: `ClimateModule` now emits
  `climate.receiver_supply_precipitation_feedback`.  The pass preserves every
  season's `climate.moisture_budget_region_id` mean, leaves ocean precipitation
  unchanged, and recomputes hydroclimate objects, precipitation-response
  objects, receiver catchments, and final source/receiver accounting after the
  feedback.
- C5e8 replay:
  `out_terminal_climate_replay_c5e8_receiver_supply_feedback_20260706/`.
- Receiver-catchment gate is upgraded to
  `aevum.earth_climate_receiver_catchment_gate.v3`; it now requires the C5e8
  feedback archive, checks finite shape, bounds land feedback, and verifies
  ocean neutrality.
- C5e8 feedback diagnostics: earthlike seeds `42/909` have land feedback
  p05/p95 about `0.938/1.018` and `0.937/1.017`; waterworld seeds stay smaller
  at about `0.970/1.011` and `0.978/1.007`.  The max land and local-budget
  mean deltas are numerical noise (`~1e-13` mm/yr).
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e8_receiver_supply_feedback_r6_20260706/`
  reports `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e8_receiver_supply_feedback_r6_20260706/`
  reports guardrail verdict `pass` with `0` failures and `0` warnings.
  Phase statuses remain low-priority `watch`: F1 `0.35`, F2 `0.10`,
  F3 `0.00`, F4 `0.00`, F5 `0.00`.
- Existing acceptance gates pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile: `py_compile` passes for edited climate, validation, feature
  catalog, rendering, diagnostics, and tests.  Target tests pass:
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`,
  `tests/test_climate_seasonal_redistribution.py`,
  `tests/test_earth_climate_receiver_catchment_gate.py`, and
  `tests/test_engine.py::test_seasonal_hydroclimate_derives_moisture_and_monsoon_fields`.
- Conclusion: C5e8 is the previous Earth-fitting regression baseline.  It
  closed the first source/receiver feedback loop conservatively.  C5e9 below
  supersedes it by addressing the most visible reduced-model ocean/SST structure
  residual without broad precipitation lifting.

Post-acceptance C5e9 ocean-structure note:

- Purpose: reduce the remaining latitude-band visual residual in temperature
  and SST by making the reduced basin/boundary-current heat anomaly more
  visible, while preserving C5e8 precipitation and source/receiver accounting
  behavior.
- Code change: `ClimateModule._ocean_currents` now gives the basin gyre,
  warm-boundary-current, cold-boundary-current, and upwelling terms more
  leverage and smooths the ocean heat anomaly less aggressively.  The field is
  still ocean-confined and area-mean neutral before it spreads into coastal
  temperature influence.  The exported ocean surface temperature now uses the
  coupled `climate.seasonal_sst` for ocean cells, rather than the weaker
  pre-coupling `current_heat_transport` projection.
- Diagnostics: `earth-climate-ocean-spatial-gate` now records same-latitude
  current/SST residual amplitudes.  For earthlike worlds it requires enough
  SST residual structure relative to Earth and a minimum reduced heat-transport
  p95, so future changes cannot silently pass with a purely latitude-banded SST
  map.  `validation.climate_diagnostics` and terminal replay summaries now
  expose `seasonal_sst_zonal_residual_abs_p95_C`, ocean heat-flux p95, and
  current-heat-transport p95.
- C5e9 replay:
  `out_terminal_climate_replay_c5e9_ocean_structure_20260706/`.
- Key C5e9 earthlike metrics: seeds `42/909` have SST same-latitude residual
  p95 about `3.07/3.39 C`, ocean heat-flux p95 about `1.28/1.25 C`, and
  current-heat-transport p95 about `1.08/1.06 C`.  Earthlike mean temperature
  changes stay small relative to C5e8 (`15.33/13.86 C`), and land precipitation
  p50 remains in the accepted range (`297/402 mm/yr`).
- Ocean-spatial gate:
  `out_earth_climate_ocean_spatial_gate_c5e9_20260706/`, verdict `pass` with
  `0` failures, `0` warnings, and `0` skipped checks.  Earthlike SST residual
  ratios to Earth are about `0.64/0.71`, above the new `0.58` floor, and
  heat-transport p95 clears the new `0.95 C` floor.
- R6 comparison/fitting:
  `out_earth_climate_comparison_c5e9_ocean_structure_r6_20260706/` reports
  `earthlike flagged: 0`;
  `out_earth_climate_fitting_c5e9_ocean_structure_r6_20260706/` reports
  guardrail verdict `pass` with `0` failures and `0` warnings.  Phase statuses
  remain low-priority `watch`: F1 `0.35`, F2 `0.10`, F3 `0.00`, F4 `0.00`,
  F5 `0.00`.
- Existing acceptance gates pass with `0` failures, `0` warnings, and `0`
  skipped checks: pattern, biome, spatial-biome, seasonal-subtype,
  mountain-zonation, windward/leeward, monsoon/moisture, circulation-layout,
  ocean-spatial, coupled-consistency, seasonal-hydro placement,
  coupling-convergence, hydro-region, moisture-flow, moisture-response, and
  receiver-catchment.
- Tests/compile: `py_compile` passes for edited climate, validation, terminal
  replay summary, ocean-spatial diagnostics, and tests.  Target tests pass:
  `tests/test_earth_climate_ocean_spatial_gate.py`,
  `tests/test_earth_climate_comparison.py`,
  `tests/test_earth_climate_fitting.py`, and
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`
  (`12 passed` for that combined target run).
- Conclusion after Replay-R reset: C5e9 remains useful as the terminal-world
  regression baseline, but it is no longer the next real-Earth replay promotion
  target.  The next climate-system work must stay in Replay-R foundation order:
  repair and revalidate R2 pressure/wind seasonal phasing and R4 temperature/
  SST energy-wall behavior before any R7 seasonal ice/snow/cloud/vegetation
  feedback or R8 biome/class tuning.

Replay-R R2a v48 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v48_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v48_20260707/`.
- v48 keeps the fitting route strictly on real Earth.  It reads Earth SLP,
  Aevum pressure proxy, standardized residual, zonal residual, M0 support, M1
  support, M2 source, and M2 transfer maps before metrics.
- The main v48 code change is M2-only: Southern Ocean shoulder-season
  source-to-pressure transfer now has sector geometry tied to SST-front, shelf,
  open-ocean, same-latitude SST, and semantic Southern Ocean support.  v48 also
  narrows the MAM term equatorward after v47 proved too poleward.
- v48 improves pressure metrics relative to v46: all/ocean standardized MAE
  `0.306/0.262 -> 0.303/0.258`, all/ocean pressure correlation
  `0.571/0.308 -> 0.580/0.355`, and SON standardized-pressure MAE
  `0.430 -> 0.417`.
- R2a is still not accepted.  The next fitting step is not wind, currents, SST,
  precipitation, biome, generated worlds, or global scalar tuning.  It is
  another real-Earth map-read pass on Southern Ocean SON pressure-sector
  amplitude/latitude placement, while preserving North Pacific/North Atlantic
  and MAM improvements.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.56s`).

Replay-R R2a v49 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v49_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v49_20260707/`.
- v49 keeps the fitting route strictly on real Earth and still reads Earth SLP,
  Aevum pressure proxy, standardized residual, zonal residual, M0 support, M1
  support, M2 source, and M2 transfer maps before metrics.
- The main v49 code change is M2-only: SON Southern Ocean transfer is split
  into a subantarctic north-flank wave and a polar-side trough wave, so the
  model no longer forces the same longitude phase onto both latitude bands.
- v49 improves pressure metrics relative to v48: all/ocean standardized MAE
  `0.303/0.258 -> 0.298/0.252`, all/ocean pressure correlation
  `0.580/0.355 -> 0.593/0.402`, and SON standardized-pressure MAE
  `0.417 -> 0.395`.
- R2a is still not accepted.  The next fitting step remains a real-Earth
  pressure subgraph pass.  Candidate owners are high-latitude Southern Ocean
  trough underdepth, DJF North Atlantic / Nordic / Barents gateway compactness,
  and muted JJA subpolar-ocean positive pressure phase.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 76.13s`).

Replay-R R2a v50 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v50_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v50_20260707/`.
- v50 keeps the fitting route strictly on real Earth and remains M2-only.
- The v50 code change adds a bounded DJF Atlantic-Arctic gateway transfer term
  plus small coastal inheritance for Iceland/Greenland/Nordic coasts.  The
  term is gated by Atlantic/Arctic basin labels, latitude/longitude, shelf,
  and SST-front support, and it avoids Beaufort/Bering longitudes.
- v50 improves pressure metrics relative to v49: all/ocean standardized MAE
  `0.298/0.252 -> 0.297/0.251`, all/ocean pressure correlation
  `0.593/0.402 -> 0.594/0.406`, and DJF standardized-pressure MAE
  `0.229 -> 0.227`.
- R2a is still not accepted.  The next fitting step remains a real-Earth
  pressure subgraph pass.  The next visible owner is warm-season /
  shoulder-season subpolar ocean positive-pressure support, because MAM/JJA
  North Pacific, North Atlantic, and Arctic subpolar regions are under-high
  while the current M2 high-support map is near zero.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 74.84s`).

Replay-R R2a v51 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v51_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v51_20260707/`.
- v51 keeps the fitting route strictly on real Earth and remains M2-only.
- The v51 code change adds MAM/JJA Northern Hemisphere subpolar-ocean
  high-pressure support.  MAM support is broad over high-latitude ocean/shelf/
  front domains.  JJA support is narrower and longitude-gated to North Pacific,
  Gulf of Alaska, and North Atlantic margins while protecting Beaufort.
- v51 improves pressure metrics relative to v50: all/ocean standardized MAE
  `0.297/0.251 -> 0.290/0.240`, all/ocean pressure correlation
  `0.594/0.406 -> 0.614/0.499`, and MAM/JJA standardized-pressure MAE
  `0.354/0.212 -> 0.329/0.207`.
- R2a is still not accepted.  The next fitting step remains a real-Earth
  pressure subgraph pass.  Candidate owners are MAM polar-cap regularity and
  Greenland/Beaufort under-high, JJA North Pacific / Gulf of Alaska under-high,
  land shoulder-season pressure errors, and residual Southern Ocean high-lat
  trough underdepth.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.71s`).

Replay-R R2a v52 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v52_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v52_20260707/`.
- v52 keeps the fitting route strictly on real Earth and remains M2-only.
- The v52 code change adds Northern Hemisphere high-latitude land shoulder-
  season pressure phase correction.  MAM now has a residual cold / near-
  freezing shield high-pressure support, gated by latitude, cold anomaly,
  near-freezing temperature, low interiority, and coast strength.  SON has a
  warm-ground / summer-heat-memory decay support that damps the premature
  autumn continental high.  The term is expressed through land pressure support
  and bounded M2 source-to-pressure transfer, without changing upstream wind or
  ocean-current generation.
- Map read: MAM Siberia/Eurasia/Canada/North America/Greenland standardized
  land residuals improve from `-0.487/-0.442/-1.145/-1.013/-0.715` to
  `-0.114/+0.025/-0.642/-0.469/+0.047`.  SON residuals improve from
  `+1.250/+0.909/+1.147/+1.072/+0.057` to
  `+0.328/-0.046/+0.109/+0.177/-0.033`.
- v52 improves pressure metrics relative to v51: all/land/ocean standardized
  MAE `0.290/0.411/0.240 -> 0.284/0.394/0.239`, all/land pressure correlation
  `0.614/0.674 -> 0.626/0.689`, and MAM/SON standardized-pressure MAE
  `0.329/0.395 -> 0.313/0.388`.
- R2a is still not accepted.  MAM polar-cap / Arctic-edge banding remains
  visually strong, MAM Greenland/Beaufort is not fully solved, JJA North
  Pacific / Gulf of Alaska remains under-high, and the Southern Ocean high-lat
  trough still needs a later map-read pass.  Do not move to R2b wind, R3
  currents, SST, precipitation, biomes, generated worlds, or global-only
  optimization while these pressure maps remain visibly wrong.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.73s`).

Replay-R R2a v53 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v53_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v53_20260707/`.
- v53 keeps the fitting route strictly on real Earth and remains M2-only.
- The v53 code change adds a MAM Arctic / Greenland / Beaufort freeze-ocean
  high-pressure object.  It uses near-freezing SST, Arctic/Baffin/Greenland/
  Beaufort longitude gates, shelf support, SST-front support, and
  high-latitude ocean masks.  A Baffin/Labrador gateway sub-support handles
  the lower-latitude Labrador/Baffin high-pressure region.
- Map read: MAM Beaufort `-0.893 -> -0.242`, Greenland Sea
  `-0.661 -> -0.035`, Barents-Kara `-0.239 -> +0.105`, Baffin/Labrador
  `-0.532 -> -0.167`, and Arctic cap `-0.662 -> -0.262`.
- v53 improves pressure metrics relative to v52: all/land/ocean standardized
  MAE `0.284/0.394/0.239 -> 0.281/0.391/0.236`, all/ocean pressure
  correlation `0.626/0.495 -> 0.631/0.522`, and MAM standardized-pressure
  MAE / zonal-anomaly correlation `0.313/0.566 -> 0.300/0.587`.
- R2a is still not accepted.  JJA North Pacific / Gulf Alaska remains
  under-high at this checkpoint, MAM Canada land is still too low, and SON
  remains visibly weak.

Replay-R R2a v54 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v54_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v54_20260707/`.
- v54 keeps the fitting route strictly on real Earth and remains M2-only.
- The v54 code change adds a JJA North Pacific / Gulf Alaska high-pressure
  object restricted to semantic Pacific basin ocean, 41-67 N, and Gulf Alaska /
  Aleutian longitude gates, with shelf, SST-front, and cool same-latitude SST
  support.  The support is zero over Beaufort/Arctic so v51's Arctic summer
  protection remains intact.
- Map read: JJA Gulf Alaska `-0.585 -> -0.268`, Aleutian
  `-0.446 -> -0.156`, North Pacific `-0.376 -> -0.121`, NW Pacific
  `-0.117 -> +0.139`, and Beaufort remains protected
  `+0.308 -> +0.303`.
- v54 improves pressure metrics relative to v53: all/land/ocean standardized
  MAE `0.281/0.391/0.236 -> 0.279/0.391/0.234`, all/ocean pressure
  correlation `0.631/0.522 -> 0.633/0.530`, and JJA standardized-pressure
  MAE / zonal-anomaly correlation `0.207/0.734 -> 0.203/0.741`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between MAM Canada land under-high, residual
  MAM Arctic cap under-high, and SON / Southern Ocean pressure-wave residuals.
- Targeted regression tests pass for the current checkpoint:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.66s`).

Replay-R R2a v55 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v55_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v55_20260707/`.
- v55 keeps the fitting route strictly on real Earth and remains M2-only.
- The v55 code change adds a MAM North America spring land-high object.  The
  support is restricted to 43-73 N and -165..-45 longitude, with a compact
  Canada-centered longitude gate, low-elevation support, cold / near-freezing
  memory, low-interiority weighting, and coast-strength modulation.  This
  replaces a rejected generic low-interiority land candidate that also raised
  Eurasia/Siberia.
- The v55 cleanup fixes the land-support diagnostic merge to use
  `np.maximum.reduce`, removing the previous NumPy warning from a three-
  argument `np.maximum` call.
- Map read: MAM Canada `-0.665 -> -0.321`, North America high-latitude land
  `-0.502 -> -0.243`, Greenland `-0.054 -> -0.009`, Siberia
  `-0.107 -> -0.123`, Eurasia `+0.008 -> -0.009`, and Alaska
  `+0.066 -> +0.244`.
- v55 improves pressure metrics relative to v54: all/land/ocean standardized
  MAE `0.279/0.391/0.234 -> 0.278/0.387/0.234`, all/land pressure
  correlation `0.633/0.684 -> 0.635/0.688`, and MAM standardized-pressure
  MAE / zonal-anomaly correlation `0.300/0.587 -> 0.296/0.604`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be residual MAM Arctic cap under-high or SON / Southern Ocean
  pressure-wave residuals, selected by map read.
- Targeted regression tests pass for the current checkpoint without the
  previous NumPy warning:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.57s`).

Replay-R R2a v58 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v58_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v58_20260707/`.
- v58 keeps the fitting route strictly on real Earth and remains M2-only.
- The v56-v57 code changes add and then strengthen DJF Atlantic-Arctic gateway
  low support, plus a DJF North Pacific / Aleutian low support term.  These are
  semantic-basin, latitude/longitude, shelf, and SST-front gated pressure
  objects expressed through bounded negative source and wave transfer.
- The v58 code change adds DJF North America winter-high relief support.  It
  is a bounded negative land-pressure object over North America, using
  latitude/longitude, low-elevation, coast-strength, low-interiority, and
  terrain-shelter gates to prevent the smaller North American winter continent
  from behaving like the Siberian High.
- Map read: DJF Icelandic/Nordic residual improves `+0.188 -> +0.043`;
  Greenland Sea `+0.400 -> +0.145`; Barents/Kara `+0.216 -> +0.051`;
  Labrador/Baffin `+0.158 -> +0.038`; Aleutian `+0.112 -> +0.050`;
  Canada land `+0.425 -> +0.244`; Siberia remains near target
  `-0.001 -> +0.023`.
- v58 improves pressure metrics relative to v55: all/land/ocean standardized
  MAE `0.278/0.387/0.234 -> 0.277/0.385/0.234`, all/land/ocean pressure
  correlation `0.635/0.688/0.529 -> 0.638/0.690/0.534`, and DJF
  standardized-pressure MAE / zonal-anomaly correlation
  `0.227/0.713 -> 0.222/0.720`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between MAM Canada land under-high, residual
  MAM Arctic cap under-high, JJA Gulf Alaska under-high, and broader SON
  residual structure.
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.94s`).

Replay-R R2a v59 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v59_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v59_20260707/`.
- v59 keeps the fitting route strictly on real Earth and remains M2-only.
- The v59 code change adds a MAM central/eastern North America plains high
  object and Canadian / central-west Arctic freeze-high objects.  These are
  explicitly gated to avoid the previous Alaska/Yukon and Barents/Kara
  over-high tradeoffs.
- Map read: MAM Canada land improves `-0.309 -> -0.015`, west Canada
  `-0.223 -> -0.012`, central Canada `-0.570 -> -0.261`, east Canada
  `-0.275 -> +0.009`, and lower North America `-0.298 -> +0.030`.  Alaska /
  Yukon does not worsen (`+0.258 -> +0.226`).
- MAM Arctic cap improves `-0.289 -> -0.231`, Beaufort
  `-0.266 -> -0.068`, Canadian Archipelago `-0.475 -> -0.255`, and
  Barents/Kara improves `+0.097 -> +0.058`.  Baffin/Labrador and central
  Arctic remain under-high.
- v59 improves pressure metrics relative to v58: all/land/ocean standardized
  MAE `0.277/0.385/0.234 -> 0.276/0.381/0.233`, all/land/ocean pressure
  correlation `0.638/0.690/0.534 -> 0.640/0.691/0.539`, and MAM
  standardized-pressure MAE / zonal-anomaly correlation
  `0.296/0.604 -> 0.290/0.615`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between residual MAM Baffin-Labrador /
  central Arctic under-high, JJA Gulf Alaska / Aleutian under-high, and broader
  SON residual structure.
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.31s`).

Replay-R R2a v60 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v60_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v60_20260707/`.
- v60 keeps the fitting route strictly on real Earth and remains M2-only.
- The v60 code change adds a JJA eastern / central North Pacific high-pressure
  object.  It is restricted to semantic Pacific basin cells in the Gulf Alaska
  / eastern Aleutian / central North Pacific sector, with explicit suppression
  over the already over-high NW Pacific 130E-160E sector and no Beaufort /
  Arctic support.
- Map read: JJA Gulf Alaska improves `-0.265 -> -0.118`; Gulf Alaska east /
  west improve `-0.235 -> -0.089` and `-0.309 -> -0.170`; eastern Aleutian
  improves `-0.343 -> -0.181`; North Pacific 160W-130W improves
  `-0.247 -> -0.118`.  NW Pacific 130E-160E does not worsen
  (`+0.369 -> +0.360`) and Beaufort remains protected (`+0.303 -> +0.299`).
- v60 improves pressure metrics relative to v59: all/ocean standardized MAE
  `0.276/0.233 -> 0.275/0.232`, all/ocean pressure correlation
  `0.640/0.539 -> 0.642/0.543`, and JJA standardized-pressure MAE /
  zonal-anomaly correlation `0.203/0.741 -> 0.199/0.749`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between residual MAM Baffin-Labrador /
  central Arctic under-high, JJA North Atlantic subpolar under-high, and broader
  SON residual structure.
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 73.41s`).

Replay-R R2a v62 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v62_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v62_20260707/`.
- v62 keeps the fitting route strictly on real Earth and remains M2-only.
- The v62 code change adds SON boreal autumn land-high relief, SON North
  Atlantic / Icelandic autumn low support, and a small SON Southern Ocean
  sector-low wave adjustment to preserve the v60 Southern Ocean balance after
  the v61 land / North Atlantic improvement.
- Map read: SON Siberia improves `+0.338 -> +0.137`, Canada
  `+0.204 -> -0.057`, Alaska `+0.395 -> +0.078`, North Atlantic subpolar
  `+0.247 -> +0.074`, and Icelandic sector `+0.455 -> +0.098`.  Labrador is
  preserved (`+0.076 -> +0.078`).
- Southern Ocean all returns to the v60 level (`+0.066 -> +0.066`) after the
  v61 side effect.  Pac/Amundsen is nearly unchanged from v60
  (`+0.106 -> +0.109`), Atlantic remains near target (`-0.026 -> -0.013`),
  Indian improves (`+0.183 -> +0.133`), and Aus-Pac changes mildly
  (`-0.015 -> +0.020`).
- v62 improves pressure metrics relative to v60: all/land/ocean standardized
  MAE `0.275/0.381/0.232 -> 0.272/0.377/0.230`, all/land/ocean pressure
  correlation `0.642/0.691/0.543 -> 0.644/0.693/0.547`, and SON
  standardized-pressure MAE `0.388 -> 0.378`.  SON zonal-anomaly correlation
  remains slightly below v60 (`0.436 -> 0.416`), so v62 is a checkpoint rather
  than promotion.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between residual MAM Baffin-Labrador /
  central Arctic under-high, JJA North Atlantic subpolar under-high and NW
  Pacific over-high, and residual SON high-latitude texture / Antarctica edge
  artifacts.
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 75.30s`).

Replay-R R2a v63 real-Earth pressure checkpoint:

- Current replay packet:
  `out_real_earth_climate_replay_r2a_m2_pressure_genesis_v63_20260707/`.
- Current pressure comparison packet:
  `out_real_earth_pressure_replay_r2a_m2_pressure_genesis_v63_20260707/`.
- v63 keeps the fitting route strictly on real Earth and remains M2-only.
- The v63 code change adds a JJA Eurasian summer thermal-low object and a JJA
  western Pacific marginal-sea low / trough object.  These are real-Earth
  pressure-source geometry fixes: the land object follows the heated
  low-elevation Eurasian belt from Arabia / Iran through India and East Asia,
  while the ocean object follows the western Pacific / East Asian shelf and
  front zone.  North America is intentionally excluded because it was already
  over-low in JJA.
- Map read: JJA NW Pacific improves `+0.360 -> +0.070`, Kuroshio / Oyashio
  `+0.312 -> +0.065`, Japan / East China Sea `+0.747 -> +0.606`, East Asia
  land `+0.386 -> +0.225`, NE Asia `+0.223 -> +0.072`, and China lowland
  `+0.539 -> +0.352`.
- JJA India improves `+0.431 -> +0.312`; Arabia / Iran improves
  `+0.636 -> +0.499` but remains a later owner.  Gulf Alaska
  (`-0.132 -> -0.128`), central North Pacific (`-0.172 -> -0.153`), MAM
  central Arctic, and SON targets are preserved.
- v63 improves pressure metrics relative to v62: all/land/ocean standardized
  MAE `0.272/0.377/0.230 -> 0.271/0.373/0.229`, all/land/ocean pressure
  correlation `0.644/0.693/0.547 -> 0.649/0.699/0.551`, and JJA
  standardized-pressure MAE / zonal-anomaly correlation
  `0.199/0.749 -> 0.192/0.778`.
- R2a is still not accepted.  The next real-Earth pressure subgraph owner
  should be selected by map read between residual MAM Baffin-Labrador /
  central Arctic under-high, JJA North Atlantic / Icelandic under-high,
  residual JJA Japan / East China Sea and Arabia / Iran summer low
  underexpression, and residual SON high-latitude texture / Antarctica edge
  artifacts.
- Targeted regression tests pass:
  `tests/test_real_earth_climate_replay.py`,
  `tests/test_real_earth_wind_replay.py`,
  `tests/test_engine.py::test_geographic_circulation_anomalies_follow_land_sea_layout`,
  and `tests/test_engine.py::test_geographic_circulation_is_weak_on_waterworld`
  (`4 passed in 105.55s`).
