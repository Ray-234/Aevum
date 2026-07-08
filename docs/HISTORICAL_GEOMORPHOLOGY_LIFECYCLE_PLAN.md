# Historical Geomorphology Lifecycle Plan

Status: active implementation plan
Owner: tectonics / terrain / archive diagnostics integration
Created: 2026-07-05
Current entry point: P174 time-continuous relief response and anti-pop guards

This document archives the plan for turning late terminal-only terrain and
bathymetry quality into time-continuous geological process output.  The trigger
for this plan was the elevation-evolution video review: the `4500 Myr` endpoint
can now look plausible, but earlier archive frames still expose broad ordinary
continental plateaus and simplified ocean basins.  The same problem exists on
land and under the ocean.

The goal is not to apply final visual polish to every historical frame.  That
would make the archive causally misleading.  The goal is to make persistent
landform and ocean-floor objects evolve through time, then let terminal polish
become a small endpoint cleanup rather than the main source of visible detail.

Related documents:

- `docs/PLATE_TECTONICS_ENGINEERING_PLAN.md` - current P-series development log.
- `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md` - source, theory, benchmark,
  and real-Earth comparison archive.
- `docs/CONTINENTAL_PHYSIOGRAPHIC_ARCHITECTURE_PLAN.md` - continental province
  architecture plan.
- `docs/EARTH_GEOMORPHOLOGY_COVERAGE.md` - required geomorphology feature
  coverage contract.

## Problem Statement

Current generated maps contain two different quality regimes:

- terminal `4500 Myr` maps receive a large late stack of terrain, planform,
  orogenic, ocean-floor, and derasterization repairs;
- earlier archive frames mostly contain the raw process-layer result and do not
  receive the same detail expression.

On land, this means young, middle-aged, and even mature continents can remain
broad, high, smooth, ordinary plateau surfaces until late terminal repairs add
province diversity.  The relevant late repairs include P167/P168-style inland
plateau diversity and surface derasterization.

Under the ocean, the same pattern appears in bathymetry.  `_regionalize_ocean_floor`
does provide process-time shelves, slopes, abyssal targets, ridges, trenches,
and transforms, but P1115-style object-backed ocean-floor hierarchy expression
is explicitly late, endpoint-oriented, and Earth-like gated.  Earlier frames can
therefore miss persistent fracture zones, hotspot tracks, island chains, guyots,
microcontinents, oceanic plateaus, and fully readable ridge-age-depth structure.

## Target Outcome

Every archive frame should carry enough causally generated geomorphology to be
interpretable as a plausible state in geological time:

- early worlds may be simpler, but not featureless;
- middle-history continents should already contain cratons, platforms, rifts,
  basins, old orogens, passive margins, and active belts when their parent
  processes exist;
- middle-history oceans should already contain ridge segments, spreading-age
  bands, transform/fracture zones, trenches, island arcs, back-arc basins,
  hotspot tracks, seamount chains, oceanic plateaus, abyssal plains, and
  parented microcontinents when their parent processes exist;
- `4500 Myr` terminal maps should still improve from final polish, but disabling
  terminal polish should leave the map structurally recognizable.

## Resolution and Time-Scale Contract

Deep-time global generation uses the `8000` to `24000` cell range by default.
At this scale the model should own the high-level geological architecture:

- ocean/continent layout and open-ocean topology;
- internal continental geographic blocks and physiographic provinces;
- long-lived geological trend belts such as volcanic belts, seismic belts,
  fault zones, sutures, rifts, and mobile belts;
- first-order crustal, sedimentary, geothermal, and resource-distribution
  tendencies;
- archive-frame evolution videos and multi-world diagnostic gates.

Selected-snapshot refinement is a separate high-resolution stage, targeting
`72000` cells or higher for a chosen frame or world.  This stage may add
sub-grid and local detail such as drainage networks, rivers, lakes, erosion
landforms, coastal and submarine landforms, micro-relief, small mountain
ranges, local depressions, small islands, reefs, and atolls.

The high-resolution stage must refine the process-time parentage produced by
the global run rather than replace it.  It may densify and detail existing
provinces, basins, watersheds, shelves, islands, and ocean-floor fabrics, but
it must not silently change plate identity, ocean-basin topology, large
continental blocks, or the causal explanation archived for the frame.

## Execution Rules

- Do not solve this by replaying terminal polish on historical frames.
- Every new visible landform or ocean-floor feature must have a parent process,
  object lineage, or explicit negative-control rationale.
- Prefer persistent objects and age-dependent response over per-frame random
  texture.
- Keep process-level parallelism for multi-seed and multi-resolution validation.
  The default worker cap remains five processes unless the user changes it.
- Preserve deterministic generation for fixed seed, world type, cell count, and
  parameter set.

## Phase P170: Historical Diagnostic Baseline

Status: completed

Purpose: quantify the temporal failure before changing generation behavior.

Implementation targets:

- Add a per-frame archive audit that scans `Engine.archive.frames`.
- Emit time series metrics for land and ocean detail.
- Run the audit on the existing six evolution-video configurations:
  waterworld x2, earthlike x2, arid x2.
- Produce a compact JSON summary and optional contact sheet overlays for
  `500`, `1500`, `2500`, `3500`, and `4500 Myr`.

Land metrics:

- `inland_detail_entropy`
- `ordinary_plateau_fraction`
- `broad_flat_inland_component_count`
- `continent_province_count_p50`
- `old_orogen_expression_fraction`
- `rift_basin_expression_fraction`
- `craton_shield_platform_split_fraction`

Ocean metrics:

- `ocean_fabric_entropy`
- `ridge_visible_fraction`
- `ridge_age_symmetry_score`
- `fracture_zone_length_fraction`
- `abyssal_plain_fraction`
- `hotspot_track_count`
- `seamount_chain_count`
- `oceanic_plateau_fraction`
- `microcontinent_fraction`
- `unparented_shoal_fraction`

Acceptance:

- The audit identifies which time windows regress to ordinary plateau or
  ordinary deep-ocean expression.
- The audit makes no generation changes.
- Tests prove deterministic output and required metric keys.

## Phase P171: Persistent Geomorphology Object Layer

Status: completed

Purpose: make land and ocean features persist through time instead of being
mostly re-derived from the current surface.

Land object kinds:

- `craton_core`
- `shield`
- `stable_platform`
- `sedimentary_basin`
- `foreland_basin`
- `rift_system`
- `old_orogen`
- `young_orogen`
- `intermontane_basin`
- `passive_margin_lowland`

Ocean object kinds:

- `ridge_segment`
- `transform_fault`
- `fracture_zone`
- `abyssal_plain`
- `trench_segment`
- `island_arc`
- `backarc_basin`
- `hotspot_track`
- `seamount_chain`
- `guyot_field`
- `oceanic_plateau`
- `microcontinent`

Required object fields:

- `id`
- `kind`
- `cells`
- `birth_myr`
- `age_myr`
- `parent_process_id`
- `parent_plate_id`
- `lineage_id`
- `activity_state`
- `relief_stage`

Acceptance:

- Objects persist across archive frames unless a geological process consumes or
  transforms them.
- Object IDs are stable enough for archive diagnostics and video overlays.
- Initial implementation can be behavior-neutral; visible elevation response is
  introduced in later phases.

## Phase P172: Inland Landform Lifecycle

Status: completed

Purpose: move the useful part of terminal inland diversity into process-time
province evolution.

Process rules:

- Collision and continental assembly create `young_orogen` plus
  `foreland_basin`.
- Young orogens decay into `old_orogen` while preserving inherited linear
  boundaries.
- Stable continental interiors split into `craton_core`, `shield`, and
  `stable_platform`.
- Rifting creates `rift_system`, `rift_basin`, and `rift_shoulder` relief.
- Passive margins create `passive_margin_lowland`, coastal plains, and
  shelf-linked sediment wedges.
- Interior basins evolve by sedimentation and subsidence, not random lowering.

Implementation targets:

- stable lifecycle IDs and age/stage fields for continental landform,
  continental province, mountain-range, and plateau objects
- `_continental_detail_province`
- `_production_continental_province_graph`
- `_apply_inland_landform_region_elevation_response`
- archive fields for historical province kind, ID, age, and relief stage
- P167/P168 responsibility reduction to endpoint guardrails only

Acceptance:

- At `1500`, `2500`, `3500`, and `4500 Myr`, large continents contain multiple
  process-backed physiographic regions.
- Ordinary high plateau area falls for process reasons rather than by random
  texture.
- Old continents can be low relief without becoming featureless.

## Phase P173: Ocean-Floor Lifecycle

Status: implemented for generated-world gate; monitor in P177 visual review

Purpose: move bathymetry detail from late P1115-style endpoint expression into
continuous ocean-floor object evolution.

Process rules:

- Ridges create young oceanic crust and bilateral age bands.
- Transform offsets create persistent fracture zones and abyssal-hill fabric.
- Old oceanic crust cools and subsides into deeper abyssal plains.
- Subduction systems create trench, island arc, and back-arc basin sequences.
- Hotspots create age-progressive seamount chains, island chains, and guyots.
- LIPs create oceanic plateaus that cool, subside, erode, and may partly drown.
- Microcontinents originate from rift, terrane, or collision lineage, never from
  unsupported shallow-ocean speckle.

Implementation targets:

- `_regionalize_ocean_floor`
- `_apply_coherent_ocean_floor_fabric`
- `_ocean_fabric_objects`
- archive fields for ocean fabric kind, ID, age, parent process, ridge distance,
  and fracture-zone lineage
- P1115 responsibility reduction to endpoint readability polish only

Acceptance:

- Pre-terminal frames show readable ridges, trenches, age-depth gradients,
  fracture zones, and parented oceanic highs when their processes exist.
- Unsupported open-ocean shallow patches are either parented or deepened.
- Ocean detail remains visible in waterworld, earthlike, and arid configurations.

## Phase P174: Time-Continuous Relief Response

Status: planned

Purpose: convert persistent objects into elevation and bathymetry through
age-aware response curves.

Land response:

- `young_orogen`: high, narrow-to-moderate, continuous relief.
- `old_orogen`: lower, wider, eroded relief and natural province boundary.
- `craton_core` / `shield`: subdued relief with resistant basement texture.
- `stable_platform`: broad low to moderate relief, often sediment covered.
- `sedimentary_basin`: low, smooth, and parented by subsidence or flexure.
- `rift_system`: axial low plus shoulder uplift.
- `passive_margin_lowland`: coastal plain and shelf-linked lowland.
- broad lowland plain response: large, low-elevation, low-relief continental
  surfaces parented by stable platforms, foreland basins, intracratonic basins,
  passive-margin lowlands, drainage/sediment accumulation, or old eroded
  orogenic margins.  Earth analogues include East European Plain, North China
  Plain, and North American interior plains.  These must not be random lowered
  patches or terminal-only polish.

Ocean response:

- `ridge_segment`: shallow, continuous, age-zero axis.
- young oceanic crust: relatively shallow.
- old oceanic crust: deeper by cooling/subsidence.
- `fracture_zone`: linear relief and offset age contrast.
- `trench_segment`: narrow, deep trough.
- `seamount_chain`: age-progressive relief with subsidence and drowning.
- `oceanic_plateau`: broad shallow high with parent lineage and finite cap.
- `microcontinent`: limited continental-crust high with parented footprint.

Acceptance:

- Object relief changes smoothly across archive frames.
- Last-frame detail no longer appears suddenly.
- Earthlike generated worlds contain broad low-elevation plain/province
  systems when supported by platform, basin, passive-margin, or sedimentary
  context; continents should not be dominated by high, flat plateaus.
- Sea/land mask guards and existing final-map quality gates remain intact.

## Phase P175: Archive and Renderer Integration

Status: planned

Purpose: make the historical process layer visible and auditable.

New or strengthened archive fields:

- `terrain.historical_province_id`
- `terrain.historical_province_kind`
- `terrain.historical_province_age_myr`
- `terrain.historical_relief_stage`
- `ocean.fabric_id`
- `ocean.fabric_kind`
- `ocean.fabric_age_myr`
- `ocean.ridge_distance_myr`
- `ocean.fracture_zone_id`
- `ocean.parent_process_id`

Rendering targets:

- historical province overlay
- ocean fabric overlay
- object-lineage overlay
- elevation/bathymetry with semantic contours
- evolution videos with matching elevation and semantic panels

Acceptance:

- A historical frame can explain why a landform or ocean feature exists.
- Renderer output aligns with archive fields and terminal maps.
- Existing elevation and bathymetry rendering remains compatible.

## Phase P176: Microbenchmark Suite

Status: planned

Purpose: prevent this from becoming another visual-only tuning loop.

Required microbenchmarks:

- `P176.historical_orogen_decay`
- `P176.rift_to_passive_margin_sequence`
- `P176.craton_shield_platform_persistence`
- `P176.foreland_basin_sedimentation`
- `P176.ridge_age_symmetric_bathymetry`
- `P176.transform_fracture_zone_persistence`
- `P176.subduction_trench_arc_backarc_sequence`
- `P176.hotspot_track_age_progression`
- `P176.oceanic_plateau_subsidence`
- `P176.microcontinent_parentage`

Acceptance:

- Each benchmark isolates one geological process.
- Each benchmark has numeric thresholds and deterministic outputs.
- Failures identify the owning process layer.

## Phase P177: Multi-World Historical Visual Gate

Status: planned

Purpose: validate the lifecycle model on generated worlds, not only fixtures.

World matrix:

- waterworld x2
- earthlike x2
- arid x2

Output per world:

- elevation at `500`, `1500`, `2500`, `3500`, and `4500 Myr`
- bathymetry at the same frames
- land province/fabric semantic maps
- ocean fabric semantic maps
- object overlays
- updated evolution video

Acceptance:

- Historical frames are not visually empty before terminal time.
- Land and ocean detail are process-parented.
- Terminal maps do not regress.
- The six-world video review no longer shows a sharp terminal quality jump.

## Phase P178: Terminal Polish Contraction

Status: planned

Purpose: make late polish a guardrail rather than the source of most quality.

Terminal responsibility after contraction:

- P167/P168: local endpoint derasterization, narrow plateau guardrails, and
  preservation of already process-backed features.
- P1115: endpoint bathymetry readability, final semantic cleanup, and
  trench/ridge line clarity.

Acceptance:

- Disabling terminal polish leaves a structurally plausible map.
- Enabling terminal polish improves readability without changing geological
  meaning.
- P170/P177 diagnostics show no major pre-terminal quality cliff.

## Execution Order

1. P170: baseline diagnostics.
2. P171: persistent object schema and behavior-neutral object persistence.
3. P172: inland lifecycle and process-time elevation response.
4. P173: ocean-floor lifecycle and process-time bathymetry response.
5. P174: smooth age-aware relief response and anti-pop guards.
6. P175: archive/render integration.
7. P176: microbenchmark promotion suite.
8. P177: six-world historical visual gate.
9. P178: terminal polish contraction.

## Current Checklist

- [x] P170 historical diagnostic baseline implemented.
  - [x] Read-only per-frame archive metrics implemented.
  - [x] Archive default fields extended for historical land detail diagnostics.
  - [x] Terrain step diagnostics now record land/ocean object kind counts.
  - [x] Deterministic unit tests and small generated-world smoke completed.
  - [x] Six-world `8000`-cell waterworld/earthlike/arid baseline completed.
- [x] P171 persistent land/ocean object layer implemented.
  - [x] Archive frame object snapshots added for existing process object collections.
  - [x] P170 reads archived object snapshots before scheduler diagnostics.
  - [x] Required cross-layer object fields normalized.
  - [x] Stable object lineage/update semantics exposed and audited.
  - [x] Object persistence checked across generated-world archive frames.
- [x] P172 inland lifecycle implemented.
  - [x] Land/province/mountain/plateau objects get stable lifecycle IDs and
    standard age/stage fields across archive frames.
  - [x] First rift overpaint correction: rift-basin diagnostics and rift-margin
    sequence generation no longer treat generic shoulders/escarpments or broad
    rift-potential fields as basin coverage.
  - [x] Inland province classification calibrated for focused earthlike
    multi-seed `8000` runs so rift-basin expression no longer overpaints mature
    continental interiors in those cases.
  - [x] Age-aware relief response moved further upstream from terminal polish.
  - [x] P170/P171 generated-world gates rerun on the six-world `8000` baseline.
- [x] P173 ocean-floor lifecycle implemented.
  - [x] Ocean fabric, margin-landform, and arc/plume-landform objects get
    stable lifecycle IDs and standard age/stage fields across archive frames.
  - [x] First bounded process-time ocean-object bathymetry response added before
    terminal P1115 polish.
  - [x] P170 ocean metrics read archived ocean object masks for ridge,
    fracture-zone, plateau, microcontinent, and parented-shoal evidence.
  - [x] P173 generated-world gate promoted beyond single-seed `420` smoke.
  - [x] P173.1 frame-level unsupported-shoal attribution/debug gate added.
  - [x] P173.2 final per-frame young open-ocean age-depth floor added with
    land-mask preservation and per-frame attribution telemetry.
  - [x] P1115 responsibility contracted to endpoint readability guardrails.
- [ ] P174 time-continuous relief response implemented.
  - [x] P174.1 Large lowland plain diagnostics added for earthlike worlds.
  - [x] P174.2 Broad lowland plains generated from platform, basin,
    passive-margin, drainage, and sedimentary context.
  - [ ] P174.3 Land and ocean relief response checked for anti-pop continuity
    across archive frames.
- [x] P174.4 Selected-snapshot `72000`+ refinement remains deferred; current
    P174 work only preserves the downstream contract and must not implement
    high-resolution refinement yet.
- [ ] P175 archive and renderer integration implemented.
- [ ] P176 microbenchmarks implemented.
- [ ] P177 six-world historical visual gate passed.
- [ ] P178 terminal polish contraction completed.

## Progress Log

2026-07-05 - Plan created

- Archived the land/ocean historical-detail problem identified from the
  elevation-evolution video review.
- Set P170 as the next implementation entry point.
- Explicitly separated process-time lifecycle generation from final terminal
  polish.

2026-07-05 - P170 initial diagnostic module implemented

- Added `aevum.diagnostics.historical_geomorphology` with read-only per-frame
  land and ocean metrics for `Engine.archive.frames`.
- Added JSON writer output as `p170_historical_geomorphology_audit.json`.
- Extended `WorldArchive.DEFAULT_KEYS` with historical land-detail fields
  required by P170: continental detail, inland geomorphology region,
  continental province code, orogeny age, old-orogen decay, and rift-margin
  stage.
- Added terrain diagnostic object kind counts for continental landforms,
  margin landforms, ocean fabric, and arc/plume landforms so future archive
  frames can report ocean-object counts without replaying terminal polish.
- Added deterministic tests in `tests/test_historical_geomorphology.py`.
- Verification:
  - `python -m pytest tests/test_historical_geomorphology.py -q` -> `2 passed`.
  - Focused archive/engine/P107-related regression subset -> `10 passed,
    163 deselected`.
  - Small generated-world P170 smoke wrote
    `out_p170_smoke_20260705/p170_historical_geomorphology_audit.json`.
- Remaining P170 work: run the planned six-world `8000`-cell baseline and
  archive the resulting time-window diagnostics before marking P170 complete.

2026-07-05 - P170 six-world `8000` baseline completed

- Added `aevum.diagnostics.p170_baseline`, a reusable process-level parallel
  runner for the planned six-world baseline.
- Ran a smoke baseline:
  `out_p170_baseline_smoke_20260705/p170_six_world_baseline_summary.json`.
- Fixed `craton_shield_platform_split_fraction` so it is a true bounded
  fraction rather than a sum of overlapping shield/platform/basin masks.
  The first full baseline directory,
  `out_p170_six_world_baseline_8000_20260705`, is superseded by the corrected
  v2 baseline below.
- Ran the corrected six-world `8000`-cell baseline with `90` requested frames
  and `5` process-level workers:
  `out_p170_six_world_baseline_8000_20260705_v2/p170_six_world_baseline_summary.json`.
- Baseline acceptance:
  - `baseline_completed`: `true`
  - `six_world_8000_baseline_completed`: `true`
  - `required_metric_keys_present`: `true`
  - `generation_behavior_changed`: `false`
  - six worlds completed with `92` usable frames each.
- Coarse ordinary-plateau and ordinary-deep-ocean flags did not trigger in the
  v2 baseline, but the target-frame metrics still expose next-phase calibration
  targets: Earth-like/arid inland entropy can drop into the `0.38-0.47` range,
  rift-basin expression is often high in continental interiors, fracture-zone
  expression is sparse in many target frames, and oceanic plateau/microcontinent
  fractions are usually near zero.  These observations feed P171-P173 rather
  than more P170 threshold tuning.
- Verification:
  - `python -m pytest tests/test_historical_geomorphology.py -q` -> `2 passed`.
  - `python -m py_compile` passed for P170 diagnostic and baseline modules.

2026-07-05 - P171 archive object snapshot foundation started

- Extended `WorldArchive.Frame` with an `objects` snapshot dictionary.
- Added `WorldArchive.DEFAULT_OBJECT_KEYS` for process-relevant objects:
  tectonic boundary objects, continental provinces, continental landforms,
  margin landforms, ocean fabric, arc/plume landforms, mountain ranges,
  plateau inventory, and rift-margin sequences.
- `WorldArchive.capture()` now deep-copies those object collections alongside
  field snapshots.  This is behavior-neutral for generation but gives future
  historical diagnostics access to object lineage evidence per frame.
- Updated P170 ocean metrics to prefer archived `terrain.ocean_fabric` and
  `terrain.arc_plume_landforms` object snapshots, falling back to scheduler
  diagnostics for older archives.
- Added a P171-focused regression in `tests/test_historical_geomorphology.py`
  that verifies object snapshots are deep-copied and usable by P170.
- Verification:
  - `python -m pytest tests/test_historical_geomorphology.py -q` -> `3 passed`.
  - `python -m pytest tests/test_engine.py::test_scheduler_refreshes_state_at_final_time tests/test_historical_geomorphology.py -q`
    -> `4 passed`.
- Remaining P171 work: normalize required object fields (`birth_myr`,
  `age_myr`, `parent_process_id`, `parent_plate_id`, `lineage_id`,
  `activity_state`, `relief_stage`) and add persistence checks across generated
  archive frames.

2026-07-05 - P171 object schema normalization and persistence audit completed

- Added `P171_REQUIRED_OBJECT_FIELDS` to `aevum.archive.world_archive`.
- `WorldArchive.capture()` now normalizes archived object snapshots without
  mutating live `WorldState.objects`.  Existing aliases such as
  `formation_myr`, `formed_myr`, `start_myr`, `mean_age_myr`, `plate_id`,
  `parent_process`, `sequence_id`, `plateau_id`, and `province_id` are copied
  into the standard P171 fields when possible.
- Missing object IDs are deterministically synthesized from collection, kind,
  cells or centroid metadata, and object index.  Every archived object records
  whether the P171 required fields are present and which fields were synthesized.
- Added `aevum.diagnostics.historical_objects` with
  `historical_object_persistence_summary()` and
  `write_historical_object_audit()`.  The audit reports required-field
  completeness, collection/kind counts, unique IDs, and recurring object IDs
  across archive frames.
- Added unit coverage for deep-copy behavior, alias normalization, P170 object
  consumption, and P171 recurring-ID diagnostics.
- Generated-world smoke:
  `out_p171_object_smoke_20260705/p171_historical_object_persistence_audit.json`
  on an `800`-cell earthlike run with `12` requested frames.  Result:
  `14` archive frames, `3186` object observations, `2347` unique object IDs,
  `134` recurring object IDs, and `0` missing required field slots.
- Verification:
  - `python -m py_compile aevum/archive/world_archive.py
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/historical_objects.py
    aevum/diagnostics/p170_baseline.py
    tests/test_historical_geomorphology.py` -> passed.
  - `python -m pytest tests/test_historical_geomorphology.py -q`
    -> `4 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py -q` -> `5 passed`.

P171 remains behavior-neutral for map generation.  Its purpose is now met:
archive frames have a normalized, auditable object layer that P172/P173 can use
as the persistence contract while moving visible land and ocean relief out of
terminal-only polish.

2026-07-05 - P172 land object lifecycle stabilization started

- Added a process-object lifecycle stabilizer in `aevum.modules.terrain`.
  It matches current objects to previous-frame objects by kind, cell overlap,
  parent-process tokens, parent continent/tectonic IDs, and centroid proximity.
- Connected the stabilizer to:
  - `terrain.continental_landforms`
  - `tectonics.continental_provinces`
  - `terrain.mountain_ranges`
  - `terrain.plateau_inventory`
- Continental province objects now expose semantic `kind` values such as
  `shield`, `platform`, `intracratonic_basin`, `foreland_basin`, `active_orogen`,
  `old_orogen`, `rift_system`, `passive_margin_lowland`, and
  `volcanic_lip_plateau` instead of relying only on the generic
  `continental_province` type.
- The stabilized objects carry lifecycle fields before archive capture:
  `birth_myr`, `age_myr`, `lineage_id`, `parent_process_id`,
  `parent_plate_id`, `activity_state`, and `relief_stage`.
- Added `tests/test_terrain_lifecycle_objects.py` for ID reuse, birth/age
  preservation, deterministic new-object IDs, and required lifecycle metadata.
- Generated-world smoke:
  `out_p172_object_lifecycle_smoke_20260705/summary_compact.json` on an
  `800`-cell earthlike run with `12` requested frames.  Result: `14` archive
  frames, `3009` object observations, `2392` unique object IDs, `249`
  recurring object IDs, and `0` missing required field slots.
- The smoke also ran P170 metrics.  It passed required metric keys, but exposed
  a P172 calibration target: rift-basin expression can still overpaint much of
  the interior in some frames.  The next P172 step should tune
  `_inland_geomorphology_state`, `_continental_detail_province`, and
  `_production_continental_province_graph` so mature continental interiors split
  into shield/platform/basin/old-orogen/rift provinces more naturally.
- Follow-up diagnostic correction: P170 now only counts rift-margin stages
  `shoulder`, `rift_basin`, and `escarpment` as rift expression.  Passive-margin
  lowland/shelf/slope/rise/abyss stages no longer inflate the rift metric.
  Re-run smoke:
  `out_p172_object_lifecycle_smoke_20260705_v2/summary_compact.json`.  The same
  `800`-cell seed retained `249` recurring object IDs and `0` missing required
  field slots.  Corrected rift-basin expression still had median `1.0` and min
  `0.454`, so the next P172 step remains a real generation/classification
  calibration rather than only a diagnostic fix.
- Verification:
  - `python -m py_compile aevum/modules/terrain.py
    aevum/archive/world_archive.py aevum/diagnostics/historical_objects.py
    aevum/diagnostics/historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py` -> passed.
  - `python -m pytest tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `7 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `8 passed`.

2026-07-05 - P172 rift overpaint attribution and first calibration

- Attributed the high `rift_basin_expression_fraction` in the P172 smoke to two
  sources:
  - P170 counted `rift_shoulder` and `passive_margin_escarpment` stages as
    `rift_basin` expression.
  - `_production_bathymetry_margin_sequence()` treated broad
    `CONT_DETAIL_RIFT_BASIN` and wide `tectonics.rift_systems` influence as
    enough to write rift-margin stages through much of the inland domain.
- Tightened P170 so `rift_basin_expression_fraction` only counts
  `RIFT_MARGIN_STAGE_RIFT_BASIN`, not shoulders, escarpments, passive lowlands,
  shelves, slopes, rises, or abyssal sequence stages.
- Tightened `_production_bathymetry_margin_sequence()`:
  - `rift_basin` sequence cells now require real rift-process support rather
    than any generic low-margin `CONT_DETAIL_RIFT_BASIN` cell.
  - `rift_shoulder` no longer accepts the broad inland-state shoulder mask
    unless it lies near a real rift or passive-margin process.
  - `passive_margin_escarpment` now requires a supported rift/passive-margin
    process, so ordinary coastal-lowland fallback cannot create rift-like
    inland stages by itself.
- Tightened `_inland_geomorphology_state()`:
  - broad moderate `tectonics.rift_potential` is no longer enough to seed
    whole-continent rift provinces;
  - rift axes now require explicit rift objects, deformation axes, or high
    rift potential in low-stability/thin crust;
  - rift shoulders must stay close to those supported seeds.
- Added regression tests:
  - P170 ignores passive-lowland, shoulder, and escarpment stages for the
    rift-basin metric and counts only the true rift-basin stage.
  - P172 rift-margin sequence does not promote an entire inland domain from a
    broad synthetic shoulder candidate.
  - P172 inland state does not turn broad moderate rift potential into an
    all-rift continent without process support.
- Generated-world smoke:
  `out_p172_inland_rift_calibration_smoke_20260705/summary_compact.json` on the
  same `800`-cell earthlike seed `172` with `12` requested frames.  Result:
  `14` archive frames, `3008` object observations, `251` recurring object IDs,
  `0` missing required field slots, `0` ordinary plateau frames, and `0`
  ordinary deep-ocean frames.  `rift_basin_expression_fraction` improved to
  median `0.316` and max `0.749` from the prior corrected smoke median `0.480`
  and max `0.857`; before P170/P172 rift correction it was median `1.0`.
- Remaining P172 calibration: this is a single-seed smoke, not the full
  multi-seed `8000` gate.  The next step should check whether the residual high
  max is a real active-rift episode or another mature-interior overpaint case.
- Verification:
  - `python -m py_compile aevum/modules/terrain.py
    aevum/archive/world_archive.py aevum/diagnostics/historical_objects.py
    aevum/diagnostics/historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py` -> passed.
  - `python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py -q` -> `9 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `10 passed`.

2026-07-05 - P172 inland lifecycle gate and focused `8000` verification

- Added `aevum.diagnostics.p172_inland_lifecycle_gate`, a reusable multi-world
  diagnostic runner for P172.  It writes per-world P170/P171 audits and a
  compact aggregate summary:
  `p172_inland_lifecycle_gate_summary.json`.
- The gate reports:
  - P170 required metric completeness and ordinary-plateau frame counts;
  - P171 object-field completeness, object counts, and recurring object IDs;
  - P172 rift-basin expression max/median, mature-frame rift max, and residual
    rift-overpaint candidate labels.
- Added CLI `--job preset:label:seed` support so focused process-level parallel
  runs can be launched through `python -m` without macOS spawn failures.
- Added `tests/test_p172_inland_lifecycle_gate.py` for compact row aggregation,
  residual-rift flagging, gate acceptance counts, and CLI job parsing.
- Low-resolution smoke:
  `out_p172_inland_lifecycle_gate_smoke_420_v2_20260705/`.
  Six worlds at `420` cells, `8` requested frames, `3` workers completed with
  `0` residual-rift worlds after filtering frames with too little inland area.
  The maximum mature rift expression was `0.479`.
- Focused `8000` verification:
  `out_p172_inland_lifecycle_gate_8000_focused_20260705/`.
  Three earthlike seeds (`42`, `172`, `909`) at `8000` cells, `12` requested
  frames, `3` workers completed with:
  - `0` residual-rift worlds;
  - max mature rift expression `0.282`;
  - all P170 required metric keys present;
  - all P171 object fields complete;
  - `0` missing required object field slots in every world.
- Per-seed focused `8000` results:
  - `earthlike_seed172`: gated/mature rift max `0.133`, recurring object IDs
    `683`, ordinary plateau frames `0`.
  - `earthlike_seed42`: gated/mature rift max `0.179`, recurring object IDs
    `1410`, ordinary plateau frames `0`.
  - `earthlike_seed909`: gated/mature rift max `0.282`, recurring object IDs
    `963`, ordinary plateau frames `0`.
- P172 is not complete yet: this focused gate does not replace the planned
  six-world `8000` gate, and age-aware relief response still needs to be moved
  further upstream from terminal polish.
- Verification:
  - `python -m py_compile aevum/diagnostics/p172_inland_lifecycle_gate.py
    aevum/modules/terrain.py aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/historical_objects.py aevum/archive/world_archive.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py` -> passed.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `13 passed`.

2026-07-05 - P172 six-world `8000` inland lifecycle gate completed

- Ran the full six-world P172 gate:
  `python -m aevum.diagnostics.p172_inland_lifecycle_gate --out
  out_p172_inland_lifecycle_gate_8000_six_world_20260705 --cells 8000
  --frames 36 --max-workers 5`.
- Output:
  `out_p172_inland_lifecycle_gate_8000_six_world_20260705/p172_inland_lifecycle_gate_summary.json`.
- Acceptance:
  - `six_world_8000_gate_completed`: `true`
  - `required_metric_keys_present`: `true`
  - `object_fields_complete`: `true`
  - `p172_diagnostic_completed`: `true`
  - `generation_behavior_changed_by_gate`: `false`
- Aggregate results:
  - `0` residual rift-overpaint worlds;
  - `0` metric-incomplete worlds;
  - `0` object-field-incomplete worlds;
  - max mature rift-basin expression `0.282`;
  - median per-world rift-basin-expression max `0.244`.
- Per-world object persistence remained healthy: each world had `0` missing
  required object field slots, and recurring object IDs ranged from `692` in a
  waterworld to `5622` in an arid world.
- Per-world ordinary fallback checks stayed clean: all six worlds reported `0`
  ordinary plateau frames and `0` ordinary deep-ocean frames.
- P172 remaining work is now narrowed to moving age-aware inland relief response
  upstream from terminal polish.  That should connect stabilized lifecycle
  objects to process-time elevation so old orogens, cratons, platforms, basins,
  rift shoulders, and passive-margin lowlands appear progressively in archive
  frames instead of being primarily endpoint repairs.

2026-07-05 - P172 age-aware inland relief response completed

- Moved the inland object-age response into the existing process-time P104F
  path, before terminal polish.  `_apply_inland_landform_region_elevation_response`
  now receives stabilized `terrain.continental_landforms` and blends a bounded
  P172 object response into its province-scale correction.
- Added `_p172_age_aware_inland_lifecycle_relief_response()` in
  `aevum.modules.terrain`.  It projects stabilized landform and continental
  province objects by `kind`, `age_myr`, `birth_myr`, `relief_stage`, `cells`,
  and `province_id`; then applies conservative response curves for active
  orogens, old orogens, shields/cratons, platforms, intracratonic/foreland
  basins, rift systems, passive-margin lowlands, and volcanic/LIP plateaus.
- The response is restricted to existing inland geomorphology candidate cells,
  smoothed once, clipped to bounded meter-scale corrections, and still passes
  through P104F's land-mask preservation and elevation caps.
- Added P172 response globals:
  - `terrain.last_p172_age_aware_inland_response_object_count`
  - `terrain.last_p172_age_aware_inland_response_area_fraction`
  - `terrain.last_p172_age_aware_inland_response_mean_abs_delta_m`
  - `terrain.last_p172_age_aware_inland_response_max_abs_delta_m`
- Extended `aevum.diagnostics.p172_inland_lifecycle_gate` so future gate
  summaries record those age-aware response metrics per world and aggregate
  triggering world count, max response area fraction, and max mean absolute
  response.
- Verification:
  - `python -m py_compile aevum/diagnostics/p172_inland_lifecycle_gate.py
    aevum/modules/terrain.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py` -> passed.
  - `python -m pytest tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py -q` -> `9 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `15 passed`.
- Generated-world validation:
  - Low-resolution behavior smoke:
    `out_p172_age_aware_inland_lifecycle_gate_smoke_420_20260705/` completed
    with `0` residual rift worlds and max mature rift expression `0.479`.
  - Full six-world `8000` behavior gate:
    `out_p172_age_aware_inland_lifecycle_gate_8000_six_world_20260705/p172_inland_lifecycle_gate_summary.json`
    completed with `0` residual rift worlds, `0` metric-incomplete worlds,
    `0` object-field-incomplete worlds, max mature rift-basin expression
    `0.282`, and all six worlds reporting `0` ordinary plateau frames and
    `0` ordinary deep-ocean frames.
  - Post-summary-field smoke:
    `out_p172_age_aware_summary_metrics_smoke_420_20260705/p172_inland_lifecycle_gate_summary.json`
    confirmed the new age-aware metrics are written.  Four non-waterworld-like
    runs triggered the response; the two waterworld runs reported zero inland
    response as expected.
- P172 is now closed for the historical lifecycle plan.  The next entry point
  is P173 ocean-floor lifecycle: persistent ridge/fracture/trench/arc/hotspot/
  plateau/microcontinent objects and process-time bathymetry response.

2026-07-05 - P173 ocean-floor lifecycle first stage started

- Stabilized generation-side ocean objects with the existing lifecycle helper:
  - `terrain.ocean_fabric`
  - `terrain.margin_landforms`
  - `terrain.arc_plume_landforms`
- Extended live object lifecycle semantics so ocean objects get meaningful
  activity and relief stages such as `ridge_axis`, `active_transform`,
  `inactive_fracture_zone`, `subduction_trench`, `island_arc`,
  `backarc_extension`, `hotspot_track`, `oceanic_plateau_subsidence`, and
  `microcontinent`.
- Added `_p173_age_aware_ocean_floor_lifecycle_response()` in
  `aevum.modules.terrain`.  `_regionalize_ocean_floor()` now invokes this
  bounded process-time response before terminal P1115 polish.  The response
  consumes stabilized ocean fabric, margin, and arc/plume objects; keeps the
  original land/ocean mask; and limits object overreach by:
  - not applying broad `age_isochron` objects as another whole-basin polish
    pass;
  - leaving broad abyssal age-depth behavior to the existing crust-age response;
  - thinning overwide ridge/fracture/trench/arc/high objects to object cores;
  - clipping response magnitude to bounded meter-scale corrections.
- Added P173 response globals:
  - `terrain.last_p173_age_aware_ocean_response_object_count`
  - `terrain.last_p173_age_aware_ocean_response_area_fraction`
  - `terrain.last_p173_age_aware_ocean_response_mean_abs_delta_m`
  - `terrain.last_p173_age_aware_ocean_response_max_abs_delta_m`
- Updated P170 ocean diagnostics so archived object masks contribute to ocean
  metrics.  This makes `fracture_zone_length_fraction`,
  `oceanic_plateau_fraction`, `microcontinent_fraction`, and
  `unparented_shoal_fraction` reflect persistent object evidence instead of
  only cell-level depth/age heuristics.
- Added tests:
  - P173 ocean objects stamp ocean-specific lifecycle stage and parent fields.
  - P173 object response raises ridge/oceanic-plateau cells, deepens
    trench/fracture cores, and preserves the sea/land mask.
  - P170 reads archived ocean-fabric fracture objects into ocean metrics.
- Verification:
  - `python -m py_compile aevum/modules/terrain.py
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/historical_objects.py aevum/archive/world_archive.py
    tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py`
    -> passed.
  - `python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py -q` -> `13 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py -q` -> `17 passed`.
- Generated-world smoke:
  `out_p173_ocean_lifecycle_smoke_420_20260705_v4/summary_compact.json`.
  Earthlike seed `173`, `420` cells, `8` requested frames completed with:
  - final ocean object counts: `terrain.ocean_fabric=16`,
    `terrain.margin_landforms=15`, `terrain.arc_plume_landforms=35`;
  - P171 required fields complete with `0` missing slots and `139` recurring
    object IDs;
  - P173 response object count `15`, response area fraction `0.598`, mean
    absolute response `220.6 m`, max response `650 m`;
  - P170 object-backed ocean metrics now include fracture-zone expression
    (`fracture_zone_length_fraction` median `0.042`, max `0.112`);
  - `unparented_shoal_fraction` max stayed low at `0.00385`.
- P173 remains in progress.  The next step should promote this from a
  single-seed `420` smoke to a reusable P173 gate and then run multi-seed /
  six-world validation before contracting P1115 terminal responsibility.

2026-07-05 - P173 reusable ocean lifecycle gate added

- Added `aevum.diagnostics.p173_ocean_lifecycle_gate`, a reusable six-world
  process-level gate for waterworld, earthlike, and arid presets.  The gate
  writes per-world P170/P171 audits, summarizes ocean lifecycle collection
  counts, and checks for:
  - missing ocean lifecycle objects;
  - missing age-aware ocean response;
  - unsupported open-ocean shoal candidates;
  - overbroad P173 bathymetry response.
- The overbroad-response diagnostic now separates broad coverage from strong
  local response.  A small-area trench/ridge correction may legitimately carry
  high meter-scale amplitude, while a broad high-amplitude ocean repaint is
  still flagged.
- Added `tests/test_p173_ocean_lifecycle_gate.py` for P173 compact-summary,
  aggregate-summary, localized high-amplitude response, and CLI job parsing.
- Verification:
  - `python -m py_compile
    aevum/diagnostics/p173_ocean_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py` -> passed.
  - `python -m pytest tests/test_p173_ocean_lifecycle_gate.py
    tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py
    -q` -> `18 passed`.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py -q` -> `22 passed`.
- Generated-world validation:
  `out_p173_ocean_lifecycle_gate_smoke_420_20260705_v2/p173_ocean_lifecycle_gate_summary.json`
  completed the six-world `420`-cell, `8`-frame smoke with:
  - `0` missing ocean lifecycle object worlds;
  - `0` missing ocean response worlds;
  - `0` unsupported shoal worlds;
  - `0` overbroad ocean response worlds;
  - max unparented shoal fraction `0.00371`;
  - max ocean response area fraction `0.799`;
  - max ocean response mean absolute delta `615 m`, from a localized
    `arid_seed101` response covering only `0.011` of the world.
- P173 remains in progress.  The next execution step is an `8000` focused or
  full six-world gate, followed by P1115 responsibility contraction if the
  high-resolution diagnostics stay clean.

2026-07-05 - P173 high-resolution focused probe exposed residual mid-history shoals

- Ran the first `8000` focused P173 gate:
  `out_p173_ocean_lifecycle_gate_8000_focused_20260705/p173_ocean_lifecycle_gate_summary.json`.
  The gate completed `waterworld_seed7`, `earthlike_seed42`, and
  `earthlike_seed909` with:
  - `0` missing ocean lifecycle object worlds;
  - `0` missing ocean response worlds;
  - `0` overbroad ocean response worlds;
  - `1` unsupported shoal world: `earthlike_seed909`;
  - max unparented shoal fraction `0.0868`, peaking around `1896 Myr`.
- Added P170 attribution fixes so parented open-ocean shoals include:
  - margin landform sources such as volcanic arcs, forearc prisms, passive
    margin wedges, delta fans, and trench objects;
  - back-arc basin objects;
  - rift-margin sequence objects;
  - a one-cell object support halo matching the process-time object influence
    band used by terrain generation.
- Added a P173 process-time unsupported open-ocean shoal cleanup in
  `_regionalize_ocean_floor()`.  It deepens far-ocean shallow patches below
  `1500 m` depth only when they lack ridge/trench/transform/protected-seaway or
  ocean lifecycle object support.  Object-backed seamount, arc, plateau,
  microcontinent, and rift-margin shoals are preserved.
- Added tests covering:
  - P170 parented-shoal attribution for margin, back-arc, and rift-margin
    object halos;
  - P173 process-time shoal cleanup preserving object-backed shoals while
    sinking unsupported far-ocean shallow cells.
- Verification:
  - `python -m py_compile aevum/modules/terrain.py
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/p173_ocean_lifecycle_gate.py
    tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py
    tests/test_p173_ocean_lifecycle_gate.py` -> passed.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_terrain_lifecycle_objects.py -q` -> `24 passed`.
- Focused `earthlike_seed909` reruns:
  - after margin/back-arc attribution:
    `out_p173_ocean_lifecycle_gate_8000_earthlike909_p170_parented_fix_20260705/`
    reduced max unparented shoal fraction to `0.0831`;
  - after process-time shoal cleanup:
    `out_p173_ocean_lifecycle_gate_8000_earthlike909_process_shoal_cleanup_20260705/`
    reduced it to `0.0767`;
  - after object halo and rift-margin support:
    `out_p173_ocean_lifecycle_gate_8000_earthlike909_rift_parented_20260705/`
    reduced it to `0.0654`.
- A diagnostic rerun of the `1896 Myr` peak frame showed the remaining
  unsupported cells are open-ocean shallow abyss/rise cells rather than
  continental or plume objects:
  - median depth about `163 m`, max depth about `1121 m`;
  - median crust age about `39 Myr`;
  - `crust.type`, `crust.domain`, and `crust.origin` all oceanic/default;
  - about `88%` is `OCEAN_MARGIN_OPEN`;
  - about `88%` is `OCEAN_DEPTH_ABYSS` and about `12%` is `OCEAN_DEPTH_RISE`;
  - object support fraction of the remaining unparented mask is `0`.
- A trial that removed broad `ridge_zone` preservation made the peak worse
  (`0.0880`) and was reverted.  This suggests the residual is not solved by
  simply removing ridge protection; it needs a narrower frame-level owner for
  shallow open-ocean abyss/rise patches.
- P173 remains in progress.  Next step: add a dedicated P173.1 frame-level
  unsupported-shoal attribution/debug gate that records cleanup candidate,
  preserve, support, and post-cleanup masks per frame, then use it to decide
  whether the residual should become a new object class, a depth-province
  correction, or a stricter process-time deepening rule.

2026-07-05 - P173.1 unsupported-shoal attribution gate implemented

- Added `aevum.diagnostics.p173_unsupported_shoal_attribution`, a read-only
  frame-level diagnostic that records:
  - cleanup candidate masks;
  - structural preserve masks for ridge/trench/fracture/restricted contexts;
  - object support masks with one-cell halo for ocean lifecycle objects;
  - semantic support masks for LIP/terrane/microcontinent contexts;
  - post-cleanup residual masks, fingerprints, component summaries, dominant
    depth province / margin type / shelf-width / rift-stage categories, and an
    `owner_hint`.
- Integrated P173.1 into `aevum.diagnostics.p173_ocean_lifecycle_gate`.
  Each world now writes `p173_unsupported_shoal_attribution.json`; the gate
  summary exposes `p173_1`, aggregate max cleanup-candidate fraction, aggregate
  max post-cleanup residual fraction, and `p1731_attribution_available`.
- Added tests:
  - `tests/test_p173_unsupported_shoal_attribution.py` validates that P173.1
    separates structural preserve, object support, and true residual open-ocean
    shallow abyss/rise cells.
  - `tests/test_p173_ocean_lifecycle_gate.py` now covers P173.1 compact and
    aggregate fields.
- Verification:
  - `python -m py_compile
    aevum/diagnostics/p173_unsupported_shoal_attribution.py
    aevum/diagnostics/p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py` -> passed.
  - `python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py
    tests/test_terrain_lifecycle_objects.py -q` -> `26 passed`.
- Six-world `420` smoke:
  `out_p1731_unsupported_shoal_attribution_smoke_420_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed with:
  - `0` unsupported shoal worlds;
  - `0` overbroad ocean response worlds;
  - `0` P173.1 attribution-missing worlds;
  - max cleanup candidate fraction `0.0529`;
  - max post-cleanup residual fraction `0.0`.
- High-resolution pressure seed:
  `out_p1731_unsupported_shoal_attribution_8000_earthlike909_20260705/`.
  The `earthlike_seed909` single-world `8000`, `12`-frame gate completed with
  P173.1 attribution available but still failed unsupported-shoal acceptance:
  - max post-cleanup residual fraction `0.0654`;
  - peak at `1896.49 Myr`;
  - cleanup candidate fraction `0.127`;
  - structural preserve fraction of candidate `0.363`;
  - object support fraction of candidate `0.246`;
  - residual fraction of candidate `0.515`;
  - owner hint `open_ocean_young_shallow_abyss_rise`;
  - dominant residual depth province `abyss` at `0.883`;
  - dominant residual margin type `open` at `0.878`;
  - residual depth median `163 m`;
  - residual crust-age median `39 Myr`;
  - largest residual components are broad shallow open-ocean young-crust patches,
    not continental, plume, or lifecycle-object-backed highs.
- P173 remains in progress.  The next implementation step should target the
  newly identified owner class: young open-ocean abyss/rise cells that remain
  far too shallow after process-time cleanup.  Candidate fixes are:
  1. tighten process-time age-depth enforcement for open-ocean abyss/rise after
     object response and smoothing;
  2. add a separate transient `young_shallow_ridge_apron` object only if these
     patches are actually ridge-proximal and line-organized;
  3. correct depth-province assignment if these are mislabeled abyss/rise cells
     that should be ridge/restricted provinces.

2026-07-05 - P173.2 final per-frame young open-ocean age-depth floor promoted

- The first P173.2 attempt inserted young open-ocean age-depth enforcement inside
  `_regionalize_ocean_floor()`.  That was rejected after the `8000`
  `earthlike_seed909` pressure seed worsened: max P173.1 post-cleanup residual
  increased from `0.0654` to `0.0793`, and cleanup-candidate area increased from
  `0.127` to `0.171`.  Root cause: `_regionalize_ocean_floor()` is called before
  several sea-level, seaway, land-payback, semantic rebuild, and terminal
  bathymetry stages, so an early depth floor feeds back into planform evolution
  instead of acting as a final per-frame cleanup.
- Reworked P173.2 as `_p1732_final_young_open_ocean_depth_floor()` in
  `aevum/modules/terrain.py`.
  - It runs after final ocean/margin semantic context is available and just
    before final drainage/export packaging.
  - It preserves the current land/ocean mask; if the mask would change, the pass
    reverts.
  - It protects ridge/trench/transform regions plus current-frame ocean fabric,
    margin landforms, arc/plume landforms, rift-margin sequences, LIP/plume
    plateaus, microcontinents, and accreted terranes.
  - If it changes bathymetry, ocean semantic context and bathymetry margin
    sequence are rebuilt before export.
- Extended P173.1 attribution to report per-frame P173.2 telemetry:
  - `terrain.last_p1732_young_open_ocean_age_depth_candidate_fraction`;
  - `terrain.last_p1732_young_open_ocean_age_depth_adjusted_fraction`;
  - before/after mean depth and land-mask-preservation flags.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_terrain_lifecycle_objects.py -q`
    -> passed.
  - `./.venv/bin/python -m pytest tests/test_p173_unsupported_shoal_attribution.py
    tests/test_p173_ocean_lifecycle_gate.py -q` -> passed.
  - `./.venv/bin/python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py
    tests/test_terrain_lifecycle_objects.py -q` -> `28 passed`.
- Low-resolution six-world smoke:
  `out_p1732_final_depth_floor_smoke_420_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed with `0` unsupported shoal worlds and `0` overbroad ocean response
  worlds.
- High-resolution pressure seed with refreshed telemetry:
  `out_p1732_final_depth_floor_telemetry_8000_earthlike909_20260705/`.
  The `earthlike_seed909` single-world `8000`, `12`-frame gate passed:
  - unsupported shoal worlds: `0`;
  - overbroad ocean response worlds: `0`;
  - max P173.1 cleanup candidate fraction: `0.0782`;
  - max P173.1 post-cleanup residual fraction: `0.00199`;
  - max P173.2 young open-ocean candidate/adjusted fraction: `0.0597`;
  - peak P173.2 adjustment occurred at `1896.49 Myr`, where residual became
    `0.0`;
  - remaining peak residual is a single `12`-cell component at `3019.76 Myr`
    with residual fraction `0.00199`, below the `0.05` gate threshold.
- High-resolution six-world gate:
  `out_p1732_final_depth_floor_8000_six_world_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed in `555.3 s` with:
  - `6` completed worlds;
  - `0` unsupported shoal worlds;
  - `0` overbroad ocean response worlds;
  - max P173.1 post-cleanup residual fraction `0.00199`;
  - `six_world_8000_gate_completed = true`;
  - `unsupported_shoal_gate_passed = true`.
- P173 remains in progress, but the specific young open-ocean shallow abyss/rise
  residual owner class is now under the gate threshold.  The next P173 step is
  contracting P1115's responsibility to endpoint readability guardrails and
  making sure terminal bathymetry polish no longer owns process-time morphology.

2026-07-05 - P173.3 P1115 responsibility contracted to endpoint guardrails

- Reworked `_p1115_final_ocean_floor_hierarchy_expression()` so the default
  terminal P1115 path no longer generates process morphology.
  - Default behavior is now endpoint guardrail mode via
    `terrain.enable_p1115_endpoint_readability_guardrails = 1`.
  - P1115 keeps final semantic unsupported-shoal cleanup and narrow trench
    bathymetry/readability.
  - P1115 no longer runs `_apply_coherent_ocean_floor_fabric()` or
    `_deepen_modern_earthlike_open_ocean_shoals()` by default; those process
    relief paths remain available only when
    `terrain.enable_p1115_endpoint_readability_guardrails = 0` for A/B checks.
  - Added telemetry:
    `terrain.last_p1115_endpoint_readability_guardrail_mode`,
    `terrain.last_p1115_process_relief_enabled`, and
    `terrain.last_p1115_process_relief_adjusted_area_fraction`.
- Extended `aevum.diagnostics.p173_ocean_lifecycle_gate` to include P1115
  guardrail telemetry and acceptance:
  `p1115_terminal_process_relief_contracted`.
- Tests:
  - P1115 targeted tests:
    `tests/test_p107_audit.py::test_p1115_final_ocean_floor_expression_defaults_to_endpoint_guardrails`,
    `tests/test_p107_audit.py::test_p1115_legacy_ocean_floor_expression_requires_explicit_opt_out`,
    and related trench/helper tests -> passed.
  - Related regression:
    `./.venv/bin/python -m pytest
    tests/test_engine.py::test_scheduler_refreshes_state_at_final_time
    tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py
    tests/test_terrain_lifecycle_objects.py
    tests/test_p107_audit.py::test_p1115_final_ocean_floor_expression_defaults_to_endpoint_guardrails
    tests/test_p107_audit.py::test_p1115_legacy_ocean_floor_expression_requires_explicit_opt_out
    tests/test_p107_audit.py::test_p11163_final_trench_objects_drive_narrow_bathymetry_and_depth_province
    -q` -> `31 passed`.
- Low-resolution six-world smoke:
  `out_p1733_p1115_guardrail_contraction_smoke_420_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed with:
  - `0` unsupported shoal worlds;
  - `0` overbroad ocean response worlds;
  - `p1115_endpoint_guardrail_mode_world_count = 2`;
  - `p1115_process_relief_enabled_world_count = 0`;
  - `p1115_terminal_process_relief_contracted = true`.
- High-resolution pressure seed:
  `out_p1733_p1115_guardrail_contraction_8000_earthlike909_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed with:
  - `0` unsupported shoal worlds;
  - max P173.1 post-cleanup residual fraction `0.00199`;
  - max P173.2 young open-ocean candidate/adjusted fraction `0.0597`;
  - `p1115_endpoint_guardrail_mode_world_count = 1`;
  - `p1115_process_relief_enabled_world_count = 0`;
  - `max_p1115_process_relief_adjusted_area_fraction = 0.0`;
  - `p1115_terminal_process_relief_contracted = true`.
- High-resolution six-world gate:
  `out_p1733_p1115_guardrail_contraction_8000_six_world_20260705/p173_ocean_lifecycle_gate_summary.json`
  completed in `568.8 s` with:
  - `6` completed worlds;
  - `0` unsupported shoal worlds;
  - `0` overbroad ocean response worlds;
  - max P173.1 post-cleanup residual fraction `0.00199`;
  - `p1115_endpoint_guardrail_mode_world_count = 2`;
  - `p1115_process_relief_enabled_world_count = 0`;
  - `max_p1115_process_relief_adjusted_area_fraction = 0.0`;
  - `six_world_8000_gate_completed = true`;
  - `unsupported_shoal_gate_passed = true`;
  - `p1115_terminal_process_relief_contracted = true`.
- P173 is now considered implemented for the generated-world diagnostic gate.
  Remaining risk moves to P177 visual review: the six-world historical videos
  still need to confirm that ocean detail remains readable without a terminal
  quality cliff.

2026-07-05 - Resolution-tier and broad lowland requirements archived for P174

- Added an explicit resolution/time-scale contract:
  - `8000` to `24000` cells remain the default deep-time global generation
    scale for ocean/continent layout, internal continental blocks, geological
    trend belts, first-order resource tendencies, archive videos, and
    multi-world diagnostic gates.
  - `72000`+ cells become a selected-snapshot refinement scale for a chosen
    frame or world, adding drainage, rivers, lakes, erosion landforms, coastal
    and submarine detail, micro-relief, small islands, reefs, and atolls.
  - High-resolution refinement must inherit process-time parent IDs,
    provinces, basins, watersheds, shelves, islands, and ocean-floor fabrics
    from the global run; it must not silently change large-scale plate or ocean
    topology.
- Promoted broad low-elevation plains to a P174 acceptance target.  The model
  should generate East European Plain, North China Plain, and North American
  interior-plains analogues from stable platforms, foreland and intracratonic
  basins, passive-margin lowlands, drainage, and sedimentary context, rather
  than from random lowering or terminal-only polish.
- Split P174 into tracked subitems: lowland diagnostics, process-parented
  lowland generation, archive-frame anti-pop continuity, and preservation of
  the downstream `72000`+ refinement contract.

2026-07-05 - `72000`+ selected-snapshot refinement explicitly deferred

- Confirmed that high-resolution single-frame refinement is not part of the
  immediate P174 implementation sequence.
- Current P174 work should focus on `8000` to `24000` deep-time global maps:
  lowland diagnostics, process-parented lowland generation, and archive-frame
  continuity.
- The `72000`+ language remains only as an interface/causality constraint so
  future refinement can inherit process-time parentage instead of replacing it.

2026-07-05 - P174.1 broad lowland plain diagnostics implemented

- Extended the read-only historical geomorphology audit with large lowland
  plain metrics:
  - `lowland_plain_fraction`;
  - `broad_lowland_plain_component_count`;
  - `largest_lowland_plain_component_fraction`;
  - `lowland_plain_parented_fraction`.
- Added `lowland_plain_deficient` frame flags and summary time windows so
  generated worlds can distinguish:
  - high, flat ordinary plateaus;
  - true broad lowland plains;
  - lowland-looking cells without platform, basin, passive-margin, old-eroded,
    or sedimentary parentage.
- Added targeted tests for:
  - deterministic key completeness with the new land metrics;
  - broad parented intracratonic lowlands;
  - unparented lowlands being flagged;
  - foreland and passive-margin lowland parentage.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py -q` -> `11 passed`.
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `18 passed`.
- Generated-world smoke:
  `out_p1741_lowland_diag_smoke_720_earthlike42_20260705/` completed with
  `required_metric_keys_present = true`.
- The smoke exposed the expected current gap: final `4500 Myr` earthlike
    seed42 at `720` cells has `lowland_plain_fraction ~= 0.013`,
    `largest_lowland_plain_component_fraction ~= 0.007`, and remains
    `lowland_plain_deficient = true`.
  - This confirms P174.2 should modify process-parented lowland generation
    rather than treating P174.1 as a visual-quality pass.

2026-07-05 - P174.2 process-parented broad lowland response started

- Added a process-parented P174 lowland-plain response in the terrain layer:
  - early/process call from P104F after age-aware inland object response;
  - late preservation call before final drainage/export, after terminal
    surface derasterization and final ocean-floor cleanup;
  - parent evidence from intracratonic/foreland/rift/passive-margin basin
    provinces, stable platforms, sediment signal, platform subsidence, old
    eroded orogens, and continental landform objects;
  - land-mask preservation guard and telemetry:
    `terrain.last_p174_lowland_plain_response_area_fraction`,
    `terrain.last_p174_lowland_plain_candidate_area_fraction`,
    `terrain.last_p174_lowland_plain_parent_area_fraction`,
    `terrain.last_p174_lowland_plain_fraction_before/after`,
    `terrain.last_p174_lowland_plain_largest_component_fraction_after`,
    `terrain.last_p174_lowland_plain_parented_fraction_after`,
    `terrain.last_p174_lowland_plain_response_mean_abs_delta_m`,
    `terrain.last_p174_lowland_plain_response_land_mask_preserved`, and
    `terrain.last_p174_lowland_plain_response_stage_code`.
- Added late semantic preservation: cells adjusted by the late P174 response
  are reclassified in `terrain.continental_detail`,
  `terrain.continental_detail_region_code`, and
  `terrain.inland_geomorphology_region_code` as basin/platform lowland rather
  than remaining highland/orogen semantics.
- Corrected the P174.1 lowland diagnostic relief basis so lowland local relief
  is measured against continental-land neighbors instead of deep-ocean
  neighbors; coast-adjacent ocean depth should not invalidate a coastal or
  platform plain.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `21 passed`.
  - Related regression:
    `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `31 passed`.
- Generated-world smoke status:
  - `720`-cell earthlike seed42 remains lowland-deficient in P170 final-frame
    diagnostics despite nonzero P174 response.  Latest telemetry shows
    `terrain.last_p174_lowland_plain_response_area_fraction ~= 0.075`,
    `terrain.last_p174_lowland_plain_parent_area_fraction ~= 0.167`, and
    land-mask preservation `true`, but P170 still reports final
    `lowland_plain_fraction ~= 0.006`.
  - Root cause from mask inspection: the adjusted cells are still mostly
    narrow/coastal or adjacent to steep continental cells at this low
    resolution.  P174.2 is therefore not complete; the next step is an
    `8000`-cell validation plus additional continuity/width tuning for
    broad inland lowland belts.

2026-07-05 - P174.2 process-parented broad lowland response validated at target scale

- Preserved the user's resolution-tier decision: `72000`+ single-frame
  refinement stays deferred.  P174 validation uses the `8000` to `24000`
  deep-time global-map tier.
- Fixed historical archive coverage for the lowland diagnostics by adding
  `sediment.thickness_m` to `WorldArchive.DEFAULT_KEYS`.  P170 can now use
  sedimentary parentage in saved frames instead of treating sediment-supported
  plains as missing-evidence cases.
- Added regression coverage:
  - archive captures `sediment.thickness_m` and P170 no longer reports it as
    missing when the live world has the field;
  - a process-parented P174 sedimentary/platform lowland is recognized by P170
    as a broad parented lowland plain.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_terrain_lifecycle_objects.py
    tests/test_historical_geomorphology.py -q` -> `22 passed`.
  - `./.venv/bin/python -m pytest tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `10 passed`.
- Generated-world validation:
  `out_p1742_8000_lowland_validation_earthlike42_20260705/` completed at
  `8000` cells and 12 requested archive frames.
  - Final `4500 Myr` earthlike seed42 now passes the broad-lowland diagnostic:
    `lowland_plain_fraction ~= 0.105`,
    `broad_lowland_plain_component_count = 6`,
    `largest_lowland_plain_component_fraction ~= 0.045`, and
    `lowland_plain_parented_fraction ~= 0.966`.
  - The same run has no missing archive fields for P170.
  - `720`-cell smoke remains useful only as a fast coarse-grid warning: at
    that resolution the final continent has too few wide interior cells for
    the broad-lowland metric to be a reliable acceptance gate.
- P174.2 is considered complete for the target global resolution tier.
  Remaining work moves to P174.3 because the same `8000` run still reports
  9 lowland-deficient archive frames, mostly before the late successful
  lowland expression.  The next repair should smooth historical continuity
  and reduce time-pop in lowland relief rather than adding `72000` detail.

2026-07-05 - P174.3 anti-pop continuity diagnostics and first memory repair
started

- Added P174 continuity statistics to the P170 historical geomorphology audit:
  - mature support frame count;
  - mature lowland deficient frame count and fraction;
  - mature lowland continuity score;
  - maximum positive/negative lowland-plain fraction step;
  - terminal lowland-plain jump and terminal-pop candidate flag;
  - analogous ocean-fabric entropy step/jump fields.
- Promoted those P174 continuity fields into the P172 gate compact summary and
  aggregate output so multi-world runs can report continuity-risk labels
  directly instead of requiring manual inspection of every P170 frame row.
- Added a process-time lowland continuity memory in the terrain layer:
  `terrain.p174_lowland_plain_continuity_memory`.
  - The memory is refreshed by selected P174 lowland response cells.
  - It decays slowly on geological timescales.
  - It is cleared outside current continental, non-active-highland,
    process-supported domains.
  - It can seed later lowland response but still passes through existing
    land-mask, highland-protection, and candidate-domain guards.
- Archived `terrain.p174_lowland_plain_continuity_memory` in
  `WorldArchive.DEFAULT_KEYS` and taught P170 to treat it as explicit P174
  process parentage for lowland plains.  This fixed the intermediate failure
  mode where a lowland was visibly present but counted as unparented because
  platform/basin semantic fields had been reclassified in that frame.
- Added regression coverage for:
  - P174 continuity summary schema in P170;
  - P172 gate compact/aggregate continuity-risk reporting;
  - archive capture of P174 lowland memory;
  - P170 accepting P174 memory as lowland parentage;
  - P174 memory preserving a supported platform/basin lowland across a later
    weaker-parentage call.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `33 passed`.
- Generated-world validation:
  - `out_p1743_memory_seed_8000_earthlike42_20260705/` completed at `8000`
    cells and 12 requested archive frames.
  - Final `4500 Myr` earthlike seed42 remains good:
    `lowland_plain_fraction ~= 0.134`,
    `broad_lowland_plain_component_count = 5`,
    `largest_lowland_plain_component_fraction ~= 0.073`, and
    `lowland_plain_parented_fraction = 1.0`.
  - Parentage is substantially improved in late/mid frames because P174 memory
    is now archived and read by P170.
  - P174.3 is not complete: the same run still reports
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
- Next required step:
  - Stop using full `8000` worlds for blind tuning.
  - Add a P174.3 microbenchmark focused on broad-lowland component width,
    largest-component continuity, and local-relief smoothing across a synthetic
    platform/basin time sequence.
  - Tune the lowland candidate-width/halo/local-relief response against that
    microbenchmark, then rerun one `8000` earthlike validation.

2026-07-06 - P174.3 lowland continuity microbenchmark added

- Added a dedicated P174.3 microbenchmark:
  `test_p174_lowland_continuity_microbenchmark_bounds_component_drift`.
  The benchmark runs one synthetic platform/intracratonic-basin sequence across
  four history frames, then checks that:
  - broad lowland area remains above the acceptance floor in every frame;
  - the largest lowland component remains broad and connected;
  - P170 recognizes every frame as parented broad lowland;
  - lowland area does not drift outward frame-by-frame from P174 memory/halo
    feedback.
- Fixed the first failure exposed by that microbenchmark:
  P174 continuity memory is now refreshed only from strong process parentage
  or already-existing memory.  Plain halo cells can still be lowered in the
  current frame, but they no longer automatically become new long-lived memory.
  This prevents lowland memory from rolling outward every frame in synthetic
  platform/basin sequences.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `34 passed`.
- Generated-world validation:
  - `out_p1743_microbench_tuned_8000_earthlike42_20260706/` completed at
    `8000` cells and 12 requested archive frames.
  - Final `4500 Myr` earthlike seed42 still passes the broad-lowland endpoint:
    `lowland_plain_fraction ~= 0.127`,
    `broad_lowland_plain_component_count = 5`,
    `largest_lowland_plain_component_fraction ~= 0.068`, and
    `lowland_plain_parented_fraction ~= 0.969`.
  - P174.3 still does not pass archive continuity:
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
- Updated diagnosis:
  - Memory parentage and memory drift are no longer the main blockers.
  - The remaining generated-world failure is component geometry: several
    mature frames have total lowland area near the threshold, but the largest
    broad component is still too small or too fragmented.
- Next required step:
  - Add a second P174.3 microbenchmark for fragmented adjacent lowland
    components and component-bridge/width continuity.
  - Tune component bridging, width-aware halo, and local-relief smoothing
    against that benchmark before running another full `8000` validation.

2026-07-06 - P174.3 fragmented-component bridge microbenchmark added

- Added a second P174.3 microbenchmark:
  `test_p174_lowland_fragment_bridge_microbenchmark_connects_adjacent_components`.
  The benchmark creates several adjacent inherited lowland fragments separated
  by narrow, high, non-orogenic platform ridges.  P174 must lower/bridge the
  separators enough that P170 sees a broad, parented lowland component instead
  of several small fragments.
- Added a constrained bridge-candidate path in the P174 lowland response:
  - bridge candidates start only from existing P174 continuity-memory lowland
    source cells;
  - they require current platform context, no active-highland semantics, bounded
    relief/elevation, and platform/sediment/subsidence/memory support;
  - bridge candidates are eligible for the current frame response but are not
    automatically refreshed into long-lived memory unless they also satisfy the
    stronger memory-refresh parentage rules.
- Tried an additional selected-component bridge after component selection, but
  the `8000` earthlike validation did not improve and slightly reduced final
  parentage, so that extra path was removed.  The retained implementation is
  the pre-selection bridge-candidate path covered by the microbenchmark.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `35 passed`.
- Generated-world validation:
  - `out_p1743_fragment_bridge_8000_earthlike42_20260706/` completed at
    `8000` cells and 12 requested archive frames.
  - Final `4500 Myr` remains acceptable:
    `lowland_plain_fraction ~= 0.125`,
    `broad_lowland_plain_component_count = 6`,
    `largest_lowland_plain_component_fraction ~= 0.062`, and
    `lowland_plain_parented_fraction ~= 0.976`.
  - P174.3 still does not pass archive continuity:
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
- Updated diagnosis:
  - The local fragment-bridge failure mode is now covered by regression tests
    and has a bounded implementation.
  - The remaining generated-world problem is upstream/process-support coverage:
    early and middle mature frames still do not produce enough candidate area
    and broad largest components before P174 runs.
- Next required step:
  - Add a microbenchmark that varies stable-platform, sediment, basin, and
    inherited-memory support strength through time, then tune the lowland
    parent/candidate generation thresholds so mature supported frames reach
    the broad-lowland area and largest-component floors before terminal polish.

2026-07-06 - P174.3 stable-platform support and regional plain-response
diagnostics added

- Kept selected-snapshot `72000`+ refinement deferred.  Current work remains
  at the global `8000-24000` process-history tier.
- Added a third P174.3 microbenchmark:
  `test_p174_stable_platform_support_microbenchmark_prevents_terminal_only_lowland`.
  It varies stable-platform, high parented-basin, inherited-memory, sediment,
  and platform-subsidence support through time.  The benchmark prevents a
  failure mode where mature lowlands appear only in a final/strong-support
  frame while earlier process-supported platform frames stay as high ordinary
  plateau.
- Added P174 archive masks for future diagnostics:
  - `terrain.p174_lowland_plain_candidate_mask`
  - `terrain.p174_lowland_plain_selected_mask`
  - `terrain.p174_lowland_plain_response_mask`
- Added bounded P174 response improvements:
  - high, process-parented basin/platform catch-up candidates;
  - stage-2 broad semantic platform parentage;
  - selected-region infill for process-supported holes;
  - lowland-source bridge and local relief capping inside the same
    process-supported plain smoothing domain.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `36 passed`.
- Generated-world validation:
  - `out_p1743_stable_platform_smoke_720_earthlike42_20260706/` completed as a
    coarse smoke only.
  - `out_p1743_stable_platform_8000_earthlike42_20260706/` completed at
    `8000` cells but did not improve P174.3 continuity:
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
  - `out_p1743_region_plain_8000_earthlike42_20260706/` improved area and
    final endpoint quality but still did not pass:
    final `lowland_plain_fraction ~= 0.141`,
    final `largest_lowland_plain_component_fraction ~= 0.081`,
    `mature_lowland_deficient_frame_count = 9 / 14`.
  - `out_p1743_lowland_bridge_8000_earthlike42_20260706/` produced a small
    additional component/area improvement but still did not pass:
    final `lowland_plain_fraction ~= 0.146`,
    final `largest_lowland_plain_component_fraction ~= 0.081`,
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
- Updated diagnosis:
  - P174 can now handle the synthetic stable-platform/high-basin support
    sequence, so the local P174 response failure modes are covered better.
  - The `8000` generated-world failure has shifted from endpoint quality to
    upstream/process-support coverage and region organization:
    early mature frames still expose too little usable lowland candidate area,
    while middle frames often have enough total lowland area but the largest
    component remains below the broad-component floor.
  - Further P174-only threshold tuning shows diminishing returns.  The next
    fix should move upstream: make the process-time continental province layer
    generate broader, connected lowland/platform/basin provinces before P174
    response, instead of asking P174 to infer region objects from fragmented
    high-relief platform cells.
- Next required step:
  - Add an upstream microbenchmark for generated continental province support:
    a large mature continent should expose contiguous platform/basin/old-orogen
    lowland provinces with sufficient width before P174 runs.
  - Tune `_continental_detail_province`, inland province graph generation, and
    age-aware relief response so the P174 candidate mask starts from coherent
    region objects.  Rerun `8000` earthlike seed42 only after that microbenchmark
    passes.

2026-07-06 - P174.3 P104F/P174 telemetry bridge and generated-world failure
attribution

- Kept selected-snapshot `72000`+ refinement deferred.  Current work remains
  at the global `8000-24000` process-history tier.
- Added P170/P172 diagnostic bridge metrics for existing frame globals:
  - `p104f_pre_p174_lowland_prep_area_fraction`
  - `p104f_pre_p174_lowland_prep_mean_lowering_m`
  - `p174_lowland_plain_response_area_fraction`
  - `p174_lowland_plain_candidate_area_fraction`
  - `p174_lowland_plain_parent_area_fraction`
  - `p174_lowland_plain_continuity_memory_area_fraction`
  - `p174_lowland_plain_continuity_parent_area_fraction`
  - `p174_lowland_plain_response_mean_abs_delta_m`
  - `p174_lowland_plain_fraction_before/after`
  - `p174_lowland_plain_largest_component_fraction_after`
  - `p174_lowland_plain_parented_fraction_after`
  - `p174_lowland_plain_response_stage_code`
- Added regression coverage:
  `test_p170_land_metrics_include_p104f_p174_frame_global_telemetry`.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_p172_inland_lifecycle_gate.py tests/test_terrain_lifecycle_objects.py
    -q` -> `31 passed`.
- Generated-world validation:
  - `out_p1743_p104f_prep_8000_earthlike42_20260706/` completed at `8000`
    cells but did not pass P174.3 continuity:
    `mature_lowland_deficient_frame_count = 9 / 14`,
    `mature_lowland_continuity_score ~= 0.357`, and
    `terminal_lowland_pop_candidate = true`.
  - `out_p1743_telemetry_bridge_2000_earthlike42_20260706/` confirmed that
    the new diagnostics appear in real generated-world JSON.  It remains a
    diagnostic smoke, not an acceptance gate.
- Updated diagnosis from the telemetry bridge:
  - P104F pre-P174 lowland prep is far too narrow in generated worlds
    (`max area fraction ~= 0.015` in the `2000` smoke).
  - P174 parent support is present but not consistently converted into enough
    low, broad, connected plains (`parent area fraction` reaches roughly
    `0.05-0.17`, while candidate/response and post-response lowland fractions
    remain below the continuity floors in most mature frames).
  - The remaining failure is therefore not a missing terminal polish step.  It
    is an upstream process-object organization problem: mature continental
    platform/basin/eroded-orogen provinces are not broad and connected enough
    before P104F/P174 response, and the response targets still leave too much
    of those regions above the P170 lowland/local-relief definition.
- Next required step:
  - Add a generated-province-support microbenchmark that inspects the
    pre-P174 state directly: broad mature continents must expose connected
    platform, sag-basin, foreland/passive-margin lowland, and eroded old-orogen
    support masks with enough area and largest-component fraction.
  - Archive or compute the P104F prep mask/eligible-support mask for that
    microbenchmark, then tune upstream province generation and P104F response.
  - Rerun the `8000` earthlike seed42 validation only after that upstream
    province-support benchmark passes.

2026-07-06 - P174.3 generated-province support benchmark and local-relief
bottleneck diagnostics

- Added the generated-province-support microbenchmark:
  `test_p174_generated_province_support_infers_broad_subsiding_basin_before_p174`.
  It starts from a mature platform-dominated continent and gives the generator
  process evidence rather than hand-authored final classes: broad sediment,
  platform subsidence, stable-platform internal block context, passive-margin
  cells, and an old-subdued-orogen object.  The benchmark requires production
  province graph support for platform, broad connected intracratonic basin,
  passive margin, and old-orogen provinces before P174 runs.
- Fixed two province-graph failure modes exposed by that benchmark:
  - strong subsiding platform basins can now replace platform background
    province cells with connected intracratonic-basin provinces;
  - explicit old-orogen process objects are protected from later basin
    fallback overpaint.
- Strengthened P104F pre-P174 lowering for strong sediment/subsidence basins:
  basin prep can lower more than the former uniform `1400 m` cap while
  preserving the land/ocean mask.
- Added an initial process-lowland local-relief smoothing pass inside P174.3,
  constrained to non-active, non-LIP, non-accreted continental lowland support
  and its immediate process-backed halo.  This is deliberately bounded after
  an overly broad first attempt lowered adjacent ordinary platform too much.
- Added P170/P172 bottleneck metrics so generated-world audits can separate
  low-elevation parented support from local-relief-qualified lowland:
  - `lowland_elevation_parented_fraction`
  - `largest_lowland_elevation_parented_component_fraction`
  - `lowland_local_relief_blocked_fraction`
- Added regression coverage in
  `test_p174_lowland_plain_metrics_track_broad_parented_lowlands` for rough
  parented lowlands whose elevation support exists but local relief blocks
  formal lowland classification.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `39 passed`.
- Generated-world smoke:
  - `out_p1743_generated_province_support_2000_earthlike42_20260706/`
    completed at `2000` cells.  It did not improve the earthlike seed42
    continuity score versus `out_p1743_telemetry_bridge_2000_earthlike42_20260706`.
  - `out_p1743_process_lowland_relief_narrow_2000_earthlike42_20260706/`
    also completed at `2000` cells and still reported
    `mature_lowland_deficient_frame_count = 12 / 14` and
    `mature_lowland_continuity_score ~= 0.143`.
- Updated diagnosis:
  - The new province benchmark covers and fixes a real upstream defect, but it
    is not the dominant failure mode for the current `earthlike_seed42` smoke.
  - In-memory frame inspection showed that many mature frames already contain
    substantial process-parented low-elevation support, while P170 still sees
    too few formal broad lowland plains.  The likely blocker is local relief
    and component continuity in/after the late P174 response rather than only
    raw lowland support area.
- Next required step:
  - Run the next validation with the new P170 bottleneck metrics and inspect
    `lowland_elevation_parented_fraction` versus `lowland_plain_fraction` and
    `lowland_local_relief_blocked_fraction` per frame.
  - If the blocked fraction is high, add a narrower microbenchmark for
    post-P174/P170 local-relief qualification and tune late terrain passes so
    parented low-elevation plains remain below the `260 m` local-relief limit
    without dragging adjacent ordinary platforms down.
  - If the blocked fraction is low, redirect P174.3 work back to component
    selection/continuity-memory persistence rather than further relief
    smoothing.

2026-07-06 - P174.3 bottleneck metrics confirm local-relief blocker; harmful
apron-floor path removed

- Ran `out_p1743_bottleneck_metrics_2000_earthlike42_20260706/` with the new
  P170 bottleneck metrics.  The `2000` smoke still failed P174.3 continuity:
  `mature_lowland_deficient_frame_count = 12 / 14` and
  `mature_lowland_continuity_score ~= 0.143`.
- The new metrics clarified the dominant failure mode:
  - `lowland_elevation_parented_fraction` median ~= `0.228` of continental
    land, max ~= `0.415`;
  - `lowland_local_relief_blocked_fraction` median ~= `0.192`, max ~= `0.356`;
  - formal `lowland_plain_fraction` median remained only ~= `0.033`.
  This confirms that many process-parented low-elevation cells exist, but the
  P170 lowland classifier rejects them because adjacent continental relief is
  still too high or fragmented.
- Added the post-P174 local-relief microbenchmark:
  `test_p174_parented_lowland_relief_microbenchmark_qualifies_rough_low_plain`.
  It verifies that P174 can turn a rough, parented low-elevation basin into a
  broad P170-recognized lowland while preserving far ordinary platform.
- Tried an explicit apron-floor implementation to keep non-core platform halo
  cells above the lowland threshold.  It passed local reasoning but made the
  real `2000` smoke worse:
  `out_p1743_apron_floor_2000_earthlike42_20260706/` reported
  `mature_lowland_deficient_frame_count = 13 / 14` and
  `mature_lowland_continuity_score ~= 0.071`.
  That implementation was removed.
- Current retained code keeps only the bounded process-lowland smoothing and
  diagnostics that do not regress the focused suite.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `40 passed`.
- Updated diagnosis:
  - P174.3 is still incomplete.
  - The next implementation should not use a broad apron floor.  It should
    instead operate on selected process-lowland components directly: reduce
    internal relief and stitch adjacent same-parent lowland components without
    reclassifying/lowering ordinary platform halos into lowlands.
- Next required step:
  - Add a component-level microbenchmark using P170's new bottleneck metrics:
    large parented low-elevation components should preserve formal lowland
    interiors and largest-component area even when their edges border higher
    ordinary platforms.
  - Tune P174 selected-component smoothing/bridge logic against that benchmark,
    then rerun the `2000` bottleneck smoke before attempting another `8000`
    validation.

2026-07-06 - P174.3 P170 lowland classifier corrected to component-internal
relief

- Rechecked the `out_p1743_bottleneck_metrics_2000_earthlike42_20260706/`
  failure with an offline frame-level comparison.  The previous formal P170
  lowland classifier used the full continental neighbor relief range, so a
  valid broad low-elevation plain could be rejected whenever its boundary
  touched a higher ordinary platform or subdued highland cell.  That made the
  diagnostic over-penalize coarse-resolution lowland margins.
- Updated P170 formal lowland recognition to use local relief measured inside
  the low-elevation lowland candidate domain.  The broader continental-neighbor
  relief check remains exposed as `lowland_local_relief_blocked_fraction`, so
  audits still show where high neighboring relief is a bottleneck instead of
  silently hiding it.
- Retained the bounded P174 process-lowland smoothing and removed the harmful
  apron-floor experiment from the active code path.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `40 passed`.
- Generated-world smoke:
  - `out_p1743_internal_relief_p170_2000_earthlike42_20260706/` completed at
    `2000` cells in `141.9 s`.
  - P174 lowland continuity improved materially versus the prior bottleneck
    run:
    - `mature_lowland_deficient_frame_count`: `12 / 14` -> `5 / 14`;
    - `mature_lowland_continuity_score`: `~0.143` -> `~0.643`;
    - `lowland_plain_fraction` median: `~0.033` -> `~0.192`;
    - final `4500 Myr` `lowland_plain_fraction ~= 0.192`,
      `broad_lowland_plain_component_count = 3`, and
      `largest_lowland_plain_component_fraction ~= 0.047`.
- Updated diagnosis:
  - P174.3 has moved from severe failure to partial success.  The terminal
    map no longer looks like the only frame with meaningful lowland plains.
  - P174.3 is still not complete: deficient frames remain around
    `380-2252 Myr`, mainly because the largest formal lowland component is
    still too small or too fragmented even when total lowland area is present.
- Next required step:
  - Tune selected-component continuity and same-parent lowland stitching for
    early and middle mature frames.
  - Rerun a `2000` smoke after that targeted tuning.  Only attempt another
    `8000` validation if `mature_lowland_deficient_frame_count` drops further
    and the final-frame lowland metrics stay stable.

2026-07-06 - P174.3 selected-component stitching implemented; generated-world
gate still limited by upstream lowland-province fragmentation

- Added a targeted P174 component-level repair for lowland plains:
  - P174's internal lowland telemetry now uses relief measured inside the
    low-elevation candidate domain, while still counting only process-parented
    lowlands for its own response metrics.
  - The late lowland response now performs a constrained same-parent stitching
    pass on already low-elevation process/P170-parented cells.  It can lower
    narrow internal ribs or one-cell same-parent separators inside a lowland
    component, but it does not let newly lowered peripheral platform cells
    become fresh diagnostic anchors.
  - The earlier over-broad version that allowed current-response platform
    cells to become anchors lowered ordinary platform too much; it was
    tightened so diagnostic anchors must have been low-elevation at P174 input.
- Added regression coverage:
  - `test_p174_lowland_component_stitching_expands_parented_elevation_component`
    verifies that a large low-elevation parented component split by internal
    high ribs becomes one P170-recognized broad lowland component while far
    ordinary platform remains high.
- Tests:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `41 passed`.
- Generated-world smoke:
  - `out_p1743_component_stitch_anchor_2000_earthlike42_20260706/` completed at
    `2000` cells.
  - The gate remained unchanged from the prior P170 classifier correction:
    `mature_lowland_deficient_frame_count = 5 / 14`,
    `mature_lowland_continuity_score ~= 0.643`, and
    `terminal_lowland_pop_candidate = false`.
  - Some non-bottleneck frame metrics improved slightly, for example final
    `4500 Myr` `lowland_plain_fraction` increased from `~0.192` to `~0.196`,
    but the same early/middle frames remain deficient.
- Updated diagnosis:
  - P174 can now handle the local failure mode where a process-parented
    low-elevation component is split by internal high ribs.
  - The remaining generated-world failures are not fixed by late local
    stitching.  Failed frames still show enough total lowland area but a
    largest formal component below the gate, while P174's own response largest
    component is often smaller than the P170 low-elevation parented component.
  - This points upstream: the lowland/platform/basin province graph and P104F
    preparation are still presenting early/middle mature lowlands as several
    separate geographic provinces or weakly connected support patches before
    P174 runs.
- Next required step:
  - Add an upstream province-continuity diagnostic/microbenchmark for
    early/middle mature continents: broad platform, sag-basin, passive-margin,
    sediment, and inherited-memory support should form at least one connected
    lowland-province candidate whose area exceeds the P170 largest-component
    floor before late P174 smoothing.
  - Tune province graph/P104F prep against that benchmark before another
    `8000` validation.  Do not keep widening the late P174 apron/stitch pass
    unless the upstream diagnostic proves that connected lowland support
    already exists at P174 input.

2026-07-06 - P174.3 upstream lowland-support diagnostics and stage-2 seeding

- Added a formal upstream support contract for P174:
  - `terrain.p104f_pre_p174_lowland_support_mask` is now archived by default.
  - P170/P172 telemetry now exposes:
    - `p104f_pre_p174_lowland_support_area_fraction`;
    - `p104f_pre_p174_lowland_support_largest_component_fraction`;
    - `p104f_pre_p174_lowland_memory_seed_area_fraction`.
- Added
  `test_p174_upstream_lowland_support_forms_connected_candidate_before_p174`.
  The benchmark constructs adjacent sag-basin/platform-corridor support and
  requires P104F to expose a connected pre-P174 lowland-support component while
  preserving far ordinary platform.
- Added `_p174_seed_upstream_lowland_support(...)` and invoked it before the
  terminal stage-2 P174 lowland response.  This matters because the P170/P172
  smoke reports the final stage-2 P174 telemetry, not only the earlier P104F
  internal stage-1 response.
- Tuned upstream support recognition so it describes the process support
  domain that P174 can lower later, not only cells that are already formal
  lowlands:
  - support can include high but process-parented basins/platform corridors up
    to the P174 catch-up elevation range;
  - active orogen, plateau, LIP, and accreted-terrane domains remain excluded;
  - ordinary platform far from basin/sediment/subsidence support is not seeded.
- P174 now treats the upstream support mask as a process parent and selection
  seed, still gated by elevation, local relief, sediment/subsidence, and basin
  evidence.
- Tests:
  - `./.venv/bin/python -m py_compile aevum/modules/terrain.py
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/p172_inland_lifecycle_gate.py
    aevum/archive/world_archive.py` -> passed.
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `42 passed`.
- Generated-world smoke:
  - `out_p1743_stage2_upstream_support_2000_earthlike42_20260706/` confirmed
    that stage-2 upstream support is now visible in audit telemetry:
    support-area median rose from `0.0` to `~0.020`, and max support-largest
    component reached `~0.168` of continental area.
  - `out_p1743_support_parent_seed_2000_earthlike42_20260706/` kept the same
    P174.3 gate result:
    `mature_lowland_deficient_frame_count = 5 / 14`,
    `mature_lowland_continuity_score ~= 0.643`, and
    `terminal_lowland_pop_candidate = false`.
- Updated diagnosis:
  - The previous missing-observability problem is fixed: P170/P172 can now show
    whether connected upstream lowland support exists at the final P174 input.
  - The current generated-world bottleneck is now narrower: some failed frames
    have connected upstream support components large enough to matter, but the
    final P174 formal lowland component remains below the P170 largest-component
    floor.  For example, the `1896 Myr` frame has upstream support-largest
    `~0.058` but final formal lowland-largest `~0.022`.
  - Therefore the next repair should measure and improve support-to-plain
    conversion inside P174, not further widen support discovery.
- Next required step:
  - Add a P174 support-conversion diagnostic/microbenchmark: a connected
    upstream support component should convert into a connected formal lowland
    component above the P170 floor without lowering far ordinary platform.
  - Tune the P174 component target/relief-capping pass against that benchmark,
    then rerun `2000` before attempting another `8000` validation.

2026-07-06 - P174.3 support-to-plain conversion and halo flattening

- Added a direct support-to-plain conversion benchmark:
  - `test_p174_support_to_plain_conversion_uses_connected_upstream_support`
    builds a connected high-elevation upstream lowland-support belt with
    basin seeds and platform corridors.
  - The test requires stage-2 P174 to convert the connected support into a
    P170-recognized broad lowland component while preserving a far ordinary
    platform above `1450 m`.
- Fixed the stage-2 overreach exposed by that benchmark:
  - `broad_semantic_platform_parent` no longer treats every stable stage-2
    platform as a lowland parent.
  - Broad semantic platform support now requires upstream support, continuity
    memory, sediment/subsidence signal, lower elevation, or a high-stability
    low-relief platform context.
- Added targeted support-component conversion in P174:
  - connected `terrain.p104f_pre_p174_lowland_support_mask` components with
    process evidence are lowered to a formal lowland target;
  - an adjacent one-cell `selected & plain_smoothing_domain` halo is flattened
    so P170's neighbor-relief test does not split an otherwise valid support
    plain at the support boundary;
  - a two-cell halo was tested and rejected because it did not reduce the
    deficient-frame count and worsened the maximum negative lowland step.
- Added telemetry:
  - `p174_support_component_response_area_fraction`;
  - `p174_support_component_response_largest_component_fraction`.
  These are exposed in P170 frame land metrics and P172 compact summaries.
- Tests:
  - `./.venv/bin/python -m py_compile aevum/modules/terrain.py
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/p172_inland_lifecycle_gate.py` -> passed.
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `43 passed`.
- Generated-world smoke:
  - Effective result:
    `out_p1743_support_halo_flatten_2000f12_earthlike42_20260706/`.
  - Same-mouth comparison against
    `out_p1743_support_parent_seed_2000_earthlike42_20260706/`:
    `mature_lowland_deficient_frame_count` improved from `5 / 14` to
    `3 / 14`;
    `mature_lowland_continuity_score` improved from `~0.643` to `~0.786`;
    `terminal_lowland_pop_candidate = false`.
  - The previous multi-frame deficient windows collapsed to three isolated
    single-frame failures at approximately `380`, `1514`, and `2252 Myr`.
- Updated diagnosis:
  - The support-to-plain conversion bottleneck is partly fixed.  Large
    connected support components, notably the prior `1896 Myr` case with
    support-largest `~0.058`, now convert into P170-visible lowland structure.
  - The remaining failures are not primarily total-lowland-area failures:
    lowland fractions are already above `~0.068`, `~0.099`, and `~0.131`.
    They fail because the largest connected formal plain remains just below
    the `0.025` component floor.
  - Remaining repair should focus on process-consistent component continuity
    for these isolated early/middle frames, not on broadening ordinary
    platform lowering.
- Next required step:
  - Before promoting to `8000`, add a narrow residual-frame attribution pass
    for the three isolated failures: measure whether the split is caused by
    P170 neighbor-relief boundaries, missing P104F support, active-domain
    exclusions, or province segmentation.
  - Only then tune the smallest owner: support halo shape if it is a relief
    boundary, upstream P104F support if the support mask is too small, or
    province graph continuity if the lowland parent is semantically split.

2026-07-06 - P174.3 residual lowland-frame attribution

- Added P170 residual attribution metrics for lowland-plain failures:
  - gap metrics:
    `lowland_plain_area_gap_fraction`,
    `lowland_plain_component_gap_fraction`,
    `lowland_plain_parentage_gap_fraction`,
    `lowland_relief_boundary_gap_fraction`,
    `lowland_support_to_plain_gap_fraction`,
    `lowland_upstream_support_gap_fraction`,
    `lowland_active_exclusion_fraction`;
  - component-context metrics:
    `lowland_near_floor_plain_component_count`,
    `lowland_near_floor_parented_component_count`;
  - cause flags and dominant code:
    area-limited, parentage-limited, relief-boundary-limited,
    support-to-plain-limited, upstream-support-limited,
    active-exclusion-limited, and component-segmentation-limited.
- Added P174 continuity summary counts for mature deficient frames by cause.
  The compact P172 world summary now carries these counts directly, so future
  smoke runs do not require manual per-frame `jq` inspection.
- Added regression coverage:
  - existing parentage fixture now asserts unparented lowlands are attributed
    to parentage limitation;
  - `test_p170_lowland_residual_attribution_marks_relief_boundary_split`
    constructs a broad parented lowland split by internal high ribs and
    requires dominant attribution code `4` (`relief_boundary`).
- Tests:
  - `./.venv/bin/python -m py_compile
    aevum/diagnostics/historical_geomorphology.py
    aevum/diagnostics/p172_inland_lifecycle_gate.py` -> passed.
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `44 passed`.
- Generated-world smoke:
  - `out_p1743_residual_attribution_2000f12_earthlike42_20260706/`
    completed with the same P174.3 continuity result as the effective
    support-halo run:
    `mature_lowland_deficient_frame_count = 3 / 14`,
    `mature_lowland_continuity_score ~= 0.786`,
    `terminal_lowland_pop_candidate = false`.
  - All three remaining mature deficient frames are now attributed to
    relief-boundary limitation:
    `mature_lowland_residual_relief_boundary_frame_count = 3`,
    while area, parentage, support-to-plain, upstream-support,
    active-exclusion, and component-segmentation counts are all `0`.
  - The isolated failures remain at approximately `380`, `1514`, and
    `2252 Myr`.  Their component gaps are small (`~0.0026-0.0059`) and their
    lowland area fractions are already above the area gate, so this is a
    boundary/relief continuity problem rather than missing lowland inventory.
- Updated diagnosis:
  - Do not widen P104F upstream support or the ordinary-platform parent rules
    for the current residual failures; the attribution does not support that.
  - The next tuning owner is P174's lowland-edge relief treatment: the formal
    lowland component is still being split by high-relief boundary cells even
    when a parented low-elevation component exists above the P170 component
    floor.
  - A previous naive two-cell support halo did not help and worsened
    `lowland_plain_max_negative_step`; the next attempt should be narrower:
    target only relief-boundary cells adjacent to a near-floor parented
    lowland component, and keep the response tied to selected/support evidence.
- Next required step:
  - Add a residual relief-boundary microbenchmark that reproduces the current
    smoke failure mode: lowland area and parentage pass, the parented
    low-elevation component is large enough, but high-relief edge cells keep
    the formal plain component below `0.025`.
  - Tune P174 against that benchmark with a component-local edge cap, then run
    the same `2000f12` smoke before any `8000` validation.

2026-07-06 - P174.3 active-edge lowland semantics and top-up telemetry

- Added P170/P174 semantics for active-orogen lowland context:
  active-orogen province cells with basin/rift/sag-basin, sediment, subsidence,
  or P174-memory evidence are treated as lowland/basin edge context rather than
  generic active highland.  True `orogen`, `plateau`, volcanic-LIP, LIP,
  accreted, and high suture contexts remain excluded from lowland response.
- Added P170 active-edge buffer handling so active-orogen basin edge cells can
  count as lowland area without making the stable lowland core fail the local
  relief or broad-component checks.
- Added P174 area-limited top-up logic and telemetry:
  - `p174_area_topup_source_area_fraction`;
  - `p174_area_topup_domain_area_fraction`;
  - `p174_area_topup_candidate_area_fraction`;
  - `p174_area_topup_response_area_fraction`;
  - `p174_area_topup_response_largest_component_fraction`.
  The top-up is constrained to existing lowland/source cells and adjacent
  process-supported apron; it is not a general platform-lowering pass.
- Added tests:
  - active-orogen basin support converts to lowland while active highland stays
    high;
  - area-limited lowland source can top up from adjacent process apron;
  - P170 active-orogen basin edges no longer split broad plains.
- Verification:
  - focused suite:
    `tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py tests/test_p173_ocean_lifecycle_gate.py tests/test_p173_unsupported_shoal_attribution.py`
    passed with `48 passed`.
  - `out_p1743_active_buffer_p170_2000f12_earthlike42_20260706/` and
    `out_p1743_area_topup_telemetry_2000f12_earthlike42_20260706/` still report
    `mature_lowland_deficient_frame_count = 5 / 14` and
    `mature_lowland_continuity_score ~= 0.643`.
  - Top-up telemetry shows the remaining failed frames have very small
    top-up domains, usually nearly equal to the existing lowland source.  This
    means P174 has little adjacent process-supported apron to convert in those
    frames.
- Negative result:
  - A trial widening of upstream support corridor discovery did not improve the
    `2000f12` gate and increased support-to-plain residual attribution, so it
    was reverted.
- Current conclusion:
  - P174.3 is still incomplete.
  - The remaining `380-2252 Myr` failures are no longer best explained by a
    missing P174 halo or area top-up threshold.  The next owner should be the
    upstream continental province/support organization: early/mid continents
    need broader connected basin/platform/passive-margin support before P174
    runs, otherwise P174 can only polish isolated lowland islands.

2026-07-06 - P174.3 P170 relative-sediment semantics checked; upstream owner
confirmed

- Aligned P170 active-orogen lowland context with the P104F/P174 terrain
  semantics by computing a relative continental sediment signal from
  `sediment.thickness_m`.  Active-orogen cells with lowland-compatible detail
  and `sediment_signal >= 0.14` can now be treated as basin/lowland edge
  context even when their absolute sediment thickness is below the frame's
  high-sediment parent cutoff.
- Added a regression fixture:
  `test_p170_active_lowland_context_accepts_relative_sediment_signal`.  It
  protects the case where an active-province edge has enough relative sediment
  evidence to avoid being counted as generic active highland, while still
  preserving the broad lowland component.
- Verification:
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py -q`
    -> `12 passed`.
  - `./.venv/bin/python -m pytest tests/test_historical_geomorphology.py
    tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py
    tests/test_p173_ocean_lifecycle_gate.py
    tests/test_p173_unsupported_shoal_attribution.py -q` -> `49 passed`.
  - `out_p1743_p170_sedsignal_2000f12_earthlike42_20260706/` completed
    successfully but did not improve the generated-world continuity gate:
    `mature_lowland_deficient_frame_count = 5 / 14`,
    `mature_lowland_continuity_score ~= 0.643`,
    `terminal_lowland_pop_candidate = false`.
  - The same frames remain deficient: approximately `380`, `765`, `1139`,
    `1514`, and `2252 Myr`.  The dominant residual is still mostly
    area-limited early/mid lowland expression, with a late active-exclusion
    case at `2252 Myr`.
- Performance note:
  - The single-world `2000`-cell, `12`-frame P172/P174 smoke took about
    `23 min` on this run.  Sampling showed the worker spending time in NumPy
    array aggregation (`ufunc.at`) rather than I/O or process deadlock.  P176
    should include a smaller P174.3 upstream-support microbenchmark so this
    slow gate is reserved for promotion checks.
- Updated diagnosis:
  - The P170 relative-sediment mismatch was a real semantic gap but not the
    blocking owner for the current generated world.
  - Do not continue widening late P174 top-up/halo logic for these frames.
    The next implementation target remains upstream continental province and
    lowland-support organization: early/mid mature continents need connected
    basin/platform/passive-margin support before P174 response, with a focused
    microbenchmark that proves the connected support exists before the slow
    `2000f12` gate is rerun.

2026-07-06 - P174.3 regional apron support organization

- Added a focused upstream-support microbenchmark:
  `test_p174_upstream_support_bridges_wide_same_basin_platform_apron`.
  It builds three mature intracratonic basin seeds separated by a broad,
  weakly subsiding platform apron.  The previous local corridor-only support
  logic leaves the largest pre-P174 support component too small
  (`~0.118` of continental area), even though the cells belong to the same
  process lowland belt.
- Implemented bounded regional apron support in both P104F and
  `_p174_seed_upstream_lowland_support`:
  - it starts only from existing strong/corridor lowland support;
  - it expands through platform cells with low-to-middle relief and weak
    sediment/subsidence or low-elevation evidence;
  - it keeps orogen, plateau, active-highland, LIP, accreted, ocean, and far
    ordinary platform cells excluded.
- Verification:
  - `./.venv/bin/python -m pytest
    tests/test_terrain_lifecycle_objects.py::test_p174_upstream_support_bridges_wide_same_basin_platform_apron -q`
    -> passed.
  - `./.venv/bin/python -m pytest tests/test_terrain_lifecycle_objects.py -q`
    -> `28 passed`.
  - focused lifecycle suite:
    `tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py tests/test_p173_ocean_lifecycle_gate.py tests/test_p173_unsupported_shoal_attribution.py`
    -> `50 passed`.
  - Generated-world smoke
    `out_p1743_regional_apron_smoke_720f12_earthlike42_20260706/`:
    `mature_lowland_deficient_frame_count = 0 / 14`,
    `mature_lowland_continuity_score = 1.0`,
    `terminal_lowland_pop_candidate = false`.
- Updated diagnosis:
  - This is a positive directional fix for the upstream owner identified by
    the `2000f12` attribution, but the `720`-cell smoke is not a promotion
    gate.
  - Next required step is to rerun the slow `2000`-cell, `12`-frame earthlike
    seed42 gate.  If the mature deficient count falls below the prior `5 / 14`
    without a large negative lowland step or terminal pop, promote to an
    `8000` validation.  If it remains `5 / 14`, inspect whether the residual
    active-exclusion cells are true active highlands or lowland apron cells
    lacking archived subsidence/province evidence.

2026-07-06 - P174.3 accreted active-margin lowland context

- Ran an in-process failed-frame attribution with `tectonics.platform_subsidence`
  temporarily archived for inspection.  The `2000f12` seed42 failures were not
  dominated by LIP cells or missing platform subsidence:
  - `platform_subsidence` was effectively zero in all failed frames;
  - low-elevation excluded cells were overwhelmingly
    `DOMAIN_ACCRETED_TERRANE`, not `DOMAIN_LIP`;
  - those cells also had basin/sediment/P174-memory evidence and were often
    tagged as active-orogen or orogen detail despite being low-elevation
    sedimentary lowland context.
- Added/updated fixtures:
  - `test_p170_accreted_active_margin_basin_context_counts_as_lowland_plain`;
  - `test_p174_active_orogen_basin_support_converts_without_lowering_highland`
    now covers accreted active-margin basin support while preserving a true
    accreted active highland.
- Implemented a bounded sedimentary lowland context in P170, P104F, stage-2
  upstream support seeding, and the P174 lowland response:
  - accreted terrane cells can participate in lowland plains only when they are
    low/mid elevation and have basin, sediment, subsidence, or P174-memory
    evidence;
  - low-elevation orogen/plateau detail can be treated as lowland context only
    under the same sedimentary evidence;
  - LIP and volcanic-LIP plateau cells remain hard exclusions;
  - true high accreted/orogenic highlands remain highland-blocked.
- Verification:
  - focused lifecycle suite:
    `tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py tests/test_p173_ocean_lifecycle_gate.py tests/test_p173_unsupported_shoal_attribution.py`
    -> `51 passed`.
  - `out_p1743_accreted_lowland_2000f12_earthlike42_20260706/`:
    `mature_lowland_deficient_frame_count = 0 / 14`,
    `mature_lowland_continuity_score = 1.0`,
    `terminal_lowland_pop_candidate = false`.
  - Residual mature lowland cause counts are all `0`.
- Remaining caution:
  - `lowland_plain_max_negative_step` is now about `-0.181`, larger than the
    previous `~ -0.134`.  The deficient-frame gate is fixed for this
    `2000f12` seed42 run, but the next `8000` validation should inspect whether
    this is acceptable process fluctuation or a visible lowland-area drop.
- Next required step:
  - Promote P174.3 to an `8000` seed42 validation, then a focused multi-world
    `8000` gate if seed42 remains clean.  If the `8000` run exposes new
    residuals, first classify them by accreted-context, support-to-plain,
    area-limited, or relief-boundary owner before further tuning.

2026-07-06 - P174.3 first `8000` validation after accreted-context repair

- Ran target-tier seed42 validation:
  `out_p1743_accreted_lowland_8000f12_earthlike42_20260706/`.
- Result:
  - `mature_lowland_deficient_frame_count = 2 / 14`;
  - `mature_lowland_continuity_score ~= 0.857`;
  - `terminal_lowland_pop_candidate = false`;
  - `lowland_plain_max_negative_step ~= -0.066`;
  - `terminal_lowland_plain_jump ~= 0.150`.
- The remaining failures are now narrow target-resolution component failures,
  not missing-lowland failures:
  - `1138.6 Myr`: `lowland_plain_fraction ~= 0.159`,
    `largest_lowland_plain_component_fraction ~= 0.0224`,
    `lowland_plain_parented_fraction = 1.0`;
  - `1513.6 Myr`: `lowland_plain_fraction ~= 0.152`,
    `largest_lowland_plain_component_fraction ~= 0.0192`,
    `lowland_plain_parented_fraction = 1.0`.
  Both frames exceed the total lowland area floor by a large margin and fail
  only because the largest formal component is just below `0.025`.
- A small component-limited extension of the existing area top-up was tested in
  `out_p1743_component_topup_8000f12_earthlike42_20260706/`.  It fired but only
  produced tiny scattered top-up responses and did not improve the `2 / 14`
  result, so that code change was removed rather than kept as unproven polish.
- Current diagnosis:
  - P174.3 has moved from early/mid frames lacking valid lowland area to
    target-resolution lowland component organization: many near-floor
    parented components exist, but the response does not yet join them into a
    broad enough process plain.
  - The next owner is not accreted-context semantics, area top-up, or active
    exclusion.  It is a real component-connectivity/support graph problem at
    `8000`: build a microbenchmark with many parented lowland components and a
    process-supported corridor, then add a narrow same-parent bridge that
    increases the largest component without lowering ordinary far platform.

2026-07-06 - P174.3 component gap-repair investigation at `8000/t_end=1700`

- Implemented and tested several bounded component gap-repair attempts in
  `_p174_process_parented_lowland_plain_response`:
  - a lower component-repair gate keyed to the P174 target fraction rather than
    the full `0.105` mature lowland area;
  - local infill between existing lowland plain and process-supported repair
    candidates;
  - a bounded largest-component adjacency growth pass;
  - a P170-style diagnostic plain reference for gap-repair component labelling;
  - diagnostic-parent support for the gap-repair-only growth domain.
- Added telemetry exported through P170:
  - `p174_component_gap_repair_infill_area_fraction`;
  - `p174_component_gap_repair_largest_growth_domain_area_fraction`;
  - `p174_component_gap_repair_largest_growth_picked_area_fraction`;
  - existing gap-repair debug now records gate, candidate, accepted,
    rejection, and response fractions.
  The latest target-resolution probe below populates the largest-growth
  telemetry fields.
- Verification:
  - focused lifecycle suite:
    `tests/test_historical_geomorphology.py tests/test_terrain_lifecycle_objects.py tests/test_p172_inland_lifecycle_gate.py tests/test_p173_ocean_lifecycle_gate.py tests/test_p173_unsupported_shoal_attribution.py`
    -> `52 passed`.
  - Short target-resolution probes all remained at
    `mature_lowland_deficient_frame_count = 2 / 10`:
    - `out_p1743_gap_repair_infill_probe_8000_t1700_earthlike42_20260706/`;
    - `out_p1743_gap_repair_largest_infill_probe_8000_t1700_earthlike42_20260706/`;
    - `out_p1743_gap_repair_component_growth_probe_8000_t1700_earthlike42_20260706/`;
    - `out_p1743_gap_repair_diagnostic_plain_probe_8000_t1700_earthlike42_20260706/`;
    - `out_p1743_gap_repair_diagnostic_parent_probe_8000_t1700_earthlike42_20260706/`.
    - `out_p1743_gap_repair_growth_telemetry_probe_8000_t1700_earthlike42_20260706/`.
- Current failed-frame evidence from the latest short probe:
  - `1080.6 Myr`: P170 lowland area is already sufficient
    (`lowland_plain_fraction ~= 0.110`) and mostly parented
    (`~= 0.985`), but largest component remains `~= 0.0124`.
    Gap repair sees the P170-style plain area (`~= 0.110`) but its candidate
    area is only `~= 0.00025` of world area.  The largest-growth domain and
    picked area are both `0.0`.
  - `1295.1 Myr`: P170 lowland area is high (`~= 0.155`) and fully parented,
    but largest component remains `~= 0.0217`.  Candidate area is still tiny
    (`~= 0.00013`), and the largest-growth domain and picked area are both
    `0.0`.
- Updated diagnosis:
  - The residual is not missing total lowland area, not parentage, and not
    relief-boundary blockage.  It is a mismatch between many already-valid
    diagnostic lowland components and the much narrower process/support
    candidate graph available to P174 gap repair.
  - Local largest-component halo growth is the wrong owner for the remaining
    `8000` failures: the largest-growth domain is empty in both failed frames.
    The next implementation step should add a separate P170-aligned lowland
    component connector that traces bridge corridors between diagnostic lowland
    components through a wider diagnostic-parent/support domain, instead of
    trying to inflate one component locally.
- Performance note:
  - Each `8000/t_end=1700/8` short probe took about `4m40s`, and a full
    `8000f12` render was previously interrupted in `_regionalize_ocean_floor`.
    Before P175-P178 broad validation, add a performance pass for repeated
    dilation/regionalization hot paths or create a cheaper deterministic P174
    frame-replay fixture.
