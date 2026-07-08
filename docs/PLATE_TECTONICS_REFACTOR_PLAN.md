# Aevum Plate Tectonics Refactor Plan

Status: active long-horizon refactor plan
Owner: tectonics / terrain / archive / diagnostics integration
Created: 2026-06-22

This document tracks the tectonics-system refactor from a usable causal proxy
toward a process-model architecture.  The current engine should remain
runnable during the refactor.  Each stage must have microbenchmarks so parameter
changes can be tuned against explicit geological behavior instead of judged
only by final rendered maps.

Earth-derived geomorphology coverage is tracked in
`docs/EARTH_GEOMORPHOLOGY_COVERAGE.md`.  The refactor is not complete until the
model can generate the required Earth feature classes through process objects
and benchmarked lifecycles, not through visual coincidence or final terrain
patches.

## Refactor Goal

Replace random or late-map corrective decisions with principle-model decisions:

- random seeds may create small initial perturbations and parameter ensembles;
- random seeds should not decide major geological events such as continent
  birth, plate breakup, subduction initiation, plume location, or seaway
  opening;
- each major tectonic feature should have a persistent object and lifecycle;
- terrain and map assets should be derived from tectonic objects and physical
  state fields, not repaired after the fact;
- every phase should be covered by microbenchmarks that isolate one geological
  process and quantify whether the process makes sense;
- Earth-like output should be evaluated by feature-class coverage: cratons,
  platforms, basins, rifts, passive margins, active margins, orogens, plateaus,
  ridges, transforms, fracture zones, trenches, island arcs, hotspot tracks,
  shelves, abyssal plains, sedimentary margins, and ice-related landforms where
  the climate/cryosphere stack is active.

## Randomness Policy

Allowed uses:

- small perturbations in otherwise deterministic initial fields;
- parameter ensembles for calibration and sensitivity analysis;
- sub-grid texture that does not create or delete major geological objects;
- deterministic tie-breaking where two physically ranked candidates are equal.

Temporary uses:

- existing compact continent seed placement;
- current plate split seed selection;
- current plume candidate selection;
- representative event sampling for archive readability.

Forbidden end-state uses:

- choosing major plume locations directly from `rng.choice`;
- choosing plate breakup locations directly from `rng.choice`;
- deciding ocean gateway / seaway opening in the final terrain layer;
- compensating continental area by arbitrary ocean-to-continent flips;
- generating Wilson-cycle stage labels from current-frame boundary kind alone.

Every temporary use should be listed in a randomness audit with:

- file and line;
- current purpose;
- replacement physical field or process;
- microbenchmark that will confirm the replacement.

## Microbenchmark Philosophy

Microbenchmarks are small deterministic worlds or synthetic field fixtures used
to tune one process.  They are not full release gates.

Each benchmark should define:

- fixture: synthetic geometry, preset, grid size, time span, and fixed
  parameters;
- varied parameters: coefficients or thresholds under tuning;
- measured fields: state fields and objects to inspect;
- acceptance metrics: numeric ranges and monotonic relationships;
- visual assets: one small diagnostic image when visual structure matters;
- failure meaning: what model assumption failed if the benchmark fails.

Recommended command shape:

```text
.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R1 --out out_bench_r1
```

The first implementation can be a Python module under
`aevum/diagnostics/tectonics_bench.py` that writes:

- `tectonics_bench_summary.json`;
- per-suite CSV tables;
- small PNG contact sheets for visual checks.

## Benchmark Fixture Library

The refactor should build a reusable fixture library before replacing major
systems.

Core fixtures:

- `single_ridge_ocean`: one ocean basin with a controlled spreading ridge.
- `ocean_continent_convergence`: old oceanic plate converging on continental
  margin.
- `continent_continent_collision`: two continental blocks converging after an
  intervening ocean closes.
- `rifted_craton_margin`: stable craton with thermal anomaly and extensional
  stress.
- `passive_margin_maturation`: rifted basin with sedimented shelf/slope/rise.
- `ridge_transform_offsets`: ridge segments offset by transform boundaries.
- `plume_under_continent`: thermal anomaly below a moving continent.
- `supercontinent_breakup`: large continent over heat accumulation and inherited
  weakness.
- `island_arc_accretion`: intra-oceanic arc colliding with a continent.
- `earth_feature_catalog`: a coverage-oriented suite spanning the required
  feature classes in `docs/EARTH_GEOMORPHOLOGY_COVERAGE.md`.

Each fixture should run at low cell counts first, for example 800 to 2500
cells, then be validated against an 8000-cell reference run only after the
microbenchmark passes.

## R0. Randomness Audit And Benchmark Harness

Status: initial harness implemented 2026-06-22

Purpose:

- enumerate all random and heuristic decisions in tectonics and terrain;
- establish the microbenchmark runner and output schema;
- prevent new hidden random major-event logic from entering the engine.

Code targets:

- `aevum/modules/tectonics.py`
- `aevum/modules/terrain.py`
- `aevum/diagnostics/tectonics_bench.py` new
- `tests/test_engine.py` or a new `tests/test_tectonics_bench.py`

Microbenchmarks:

- `R0.randomness_inventory`
  - fixture: static source scan plus optional runtime event trace;
  - metrics: count of RNG calls by category, count of forbidden unclassified
    calls, list of temporary calls;
  - acceptance: zero unclassified RNG calls; no forbidden RNG calls in newly
    touched major-event code.
- `R0.benchmark_determinism`
  - fixture: same benchmark suite run twice with same parameters;
  - metrics: identical JSON summaries and stable object ids;
  - acceptance: byte-identical summaries for deterministic fixtures.
- `R0.parameter_ledger_smoke`
  - fixture: one small synthetic suite;
  - metrics: all varied coefficients appear in output metadata;
  - acceptance: benchmark result can be traced to parameter values.

Exit criteria:

- randomness audit document exists;
- R0 benchmark command writes JSON/CSV outputs; later visual benchmark suites
  write PNG contact sheets when map structure matters;
- CI-level smoke benchmark runs in under 30 seconds.

Initial implementation:

- Added `aevum.diagnostics.tectonics_bench`.
- Added R0 static randomness inventory and deterministic summary digest.
- Added `/Users/rayw/Projects/aevum/docs/TECTONICS_RANDOMNESS_AUDIT.md`.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R0 --out out_bench_r0_20260622`
- Latest R0 counts:
  - total RNG calls `32`;
  - allowed `9`;
  - temporary `14`;
  - forbidden end-state debt `9`;
  - unclassified `0`.
- Test:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> `2 passed in 5.44s`.

## R1. Deep Interior And Tectonic Potential Fields

Status: initial fields and plume-potential consumer implemented 2026-06-22

Purpose:

- introduce state fields that can drive plume, rift, and subduction tendencies;
- stop using random location selection for mantle plumes and rift initiation.

Process model:

- mantle heat anomaly evolves by slow diffusion and secular cooling;
- subduction zones create downwelling potential;
- large continents increase insulation and heat accumulation underneath;
- inherited weak zones focus extension;
- plume potential is a local maximum of heat anomaly and upwelling tendency.

Fields:

- `mantle.heat_anomaly`
- `mantle.upwelling_potential`
- `mantle.downwelling_potential`
- `lithosphere.thermal_thickness`
- `tectonics.rift_potential`
- `tectonics.plume_potential`

Microbenchmarks:

- `R1.heat_diffusion_stability`
  - fixture: synthetic hot spot and cold slab anomaly on a sphere;
  - metrics: heat anomaly remains finite, decays smoothly, conserves expected
    broad integral within tolerance;
  - acceptance: no checkerboard modes, no negative runaway, monotonic smoothing.
- `R1.continental_insulation`
  - fixture: one large continent over neutral mantle;
  - metrics: heat anomaly grows under continent interior faster than under
    adjacent ocean;
  - acceptance: interior-to-ocean heat anomaly contrast exceeds target range
    without producing global saturation.
- `R1.slab_downwelling`
  - fixture: active trench line beside old oceanic crust;
  - metrics: downwelling potential localizes along trench and decays away;
  - acceptance: p90 trench-zone potential > p90 far-field potential by a tuned
    factor.
- `R1.plume_trigger`
  - fixture: moving continent over rising heat anomaly;
  - metrics: plume candidates ranked by potential, not random choice;
  - acceptance: plume appears at local potential maximum and is reproducible
    across runs.

Exit criteria:

- plume selection uses `plume_potential`; done for `tectonics.plumes`;
- rift potential exists even if not yet fully consumed by plate breakup; done;
- benchmark summary records tunable thermal parameters; done.

Initial implementation:

- `InteriorModule` now produces `mantle.heat_anomaly`,
  `mantle.upwelling_potential`, `mantle.downwelling_potential`,
  `lithosphere.thermal_thickness`, `tectonics.rift_potential`, and
  `tectonics.plume_potential`.
- The proxy model keeps the global box thermal budget, then projects low-order
  heat anomalies to surface cells: continents insulate, ridges focus
  upwelling, trenches and old slabs focus downwelling, and lithospheric
  thickness responds to crust age, crust stability, and heat anomaly.
- `TectonicsModule._plume_activity` no longer uses `rng.choice`; plume heads are
  selected from deterministic local maxima of `tectonics.plume_potential` with
  a minimum angular spacing.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R1 --out out_bench_r1_20260622`
- Latest R1 metrics:
  - `R1.heat_diffusion_stability`: roughness `0.06947 -> 0.05068`;
  - `R1.continental_insulation`: interior-ocean heat contrast `0.02846`;
  - `R1.slab_downwelling`: trench/far downwelling p90 ratio `5.8577`;
  - `R1.plume_trigger`: selected plume distance from imposed thermal maximum
    `2.9406 deg`.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py tests/test_engine.py::test_interior_heat_budget_terms_are_physical tests/test_engine.py::test_interior_r1_potential_fields_are_bounded tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> `5 passed in 11.40s`.
- Full suite:
  `.venv/bin/python -m pytest -q`
  -> `41 passed in 387.06s`.

## R2. Plate Kinematics From Torque Proxies

Status: initial torque proxy and random refresh replacement implemented
2026-06-23

Purpose:

- replace random Euler-pole refresh with a force/torque proxy;
- make plate speed and rotation respond to ridge, slab, collision, and drag
  terms.

Process model:

- slab pull from trench length, oceanic age, and negative buoyancy proxy;
- ridge push from young ridge elevation and plate age gradient;
- collision resistance from continental convergence and crustal thickness;
- basal drag from mantle flow / plate velocity mismatch;
- transform friction from shear boundary length.

Microbenchmarks:

- `R2.single_slab_pull`
  - fixture: one oceanic plate with an old leading edge at a trench;
  - metrics: plate velocity points toward trench; speed increases with oceanic
    age and trench length;
  - acceptance: monotonic speed response over parameter sweep.
- `R2.ridge_push_symmetry`
  - fixture: symmetric ridge between two oceanic plates;
  - metrics: plates diverge with equal and opposite velocity components;
  - acceptance: net angular momentum near zero within tolerance.
- `R2.collision_locking`
  - fixture: two continental plates converging;
  - metrics: convergence speed decreases as collision resistance rises;
  - acceptance: collision does not continue at oceanic-subduction speed after
    thick crust contact.
- `R2.no_random_reorg_jitter`
  - fixture: repeated reorganization interval with identical boundary forces;
  - metrics: pole changes are explained by torque delta, not random jitter;
  - acceptance: no direct random pole refresh for active plates.

Exit criteria:

- plate motion update can run from torque terms; done for reorganization
  refresh;
- random pole jitter is removed or behind a disabled fallback; done for
  `_refresh_plate_motions`; `_split_large_plates` is now handled by the P20
  topology/rift/boundary score path;
- plate velocity diagnostics expose force component contributions; done via
  `r2_force_components` on plate objects and R2 benchmark summaries.

Initial implementation:

- Added `TectonicsModule._torque_proxy_plate_motions`.
- Reorganization refresh now computes deterministic plate pole/rate updates
  from slab pull, ridge push, collision locking, basal drag, and transform
  friction.  Motion and rate memory keep finite plate reorganization gradual
  instead of redrawing Euler poles.
- Collision zones now act as resistance/locking rather than an active force
  that pushes continents apart.
- `TectonicsModule._refresh_plate_motions` no longer uses
  `random_unit_vectors` or `rng.normal`.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R2 --out out_bench_r2_20260623`
- Latest R2 metrics:
  - `R2.single_slab_pull`: velocity dot toward trench `0.99981`; old-ocean
    rate `0.00491` > young-ocean rate `0.00473`;
  - `R2.ridge_push_symmetry`: west/east rates `0.00482 / 0.00482`, rate
    asymmetry `0.00159`, net torque ratio `0.05094`;
  - `R2.collision_locking`: thick collision rate `0.00300` <
    subduction rate `0.00488`;
  - `R2.no_random_reorg_jitter`: repeated summaries are digest-identical.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py tests/test_engine.py::test_r2_plate_motion_refresh_records_force_components tests/test_engine.py::test_geography_primitives_handle_waterworld_and_arid_layouts tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> `6 passed in 13.68s`.
- Full suite:
  `.venv/bin/python -m pytest -q`
  -> `43 passed in 390.11s`.

## R3. Persistent Boundary Resolver

Status: initial persistent boundary metadata implemented 2026-06-23

Purpose:

- turn boundary cells into persistent boundary objects with geometry, polarity,
  age, adjacent plates, and lifecycle;
- make transform, subduction, ridge, collision, passive margin, and rift
  boundaries coherent over time.

Process model:

- boundary type is resolved from relative motion, crust type, inherited
  weakness, oceanic age, and prior boundary identity;
- subduction polarity chooses the older/denser or weaker descending plate;
- transform boundaries connect ridge offsets and oblique contacts;
- collision zones broaden into orogens instead of remaining one-cell sutures.

Microbenchmarks:

- `R3.subduction_polarity_old_ocean`
  - fixture: old oceanic plate vs young oceanic plate convergence;
  - metrics: old plate selected as descending side;
  - acceptance: polarity flips when age contrast is reversed.
- `R3.ocean_continent_margin`
  - fixture: ocean-continent convergence;
  - metrics: trench offshore, active margin on continent side, arc behind
    trench;
  - acceptance: trench is not painted through continental interior.
- `R3.transform_from_offset_ridge`
  - fixture: two offset spreading segments;
  - metrics: transform object links ridge tips and has shear-dominant motion;
  - acceptance: transform length and orientation match imposed offset.
- `R3.boundary_persistence`
  - fixture: stable ridge/trench through several frames;
  - metrics: object id persistence, geometry overlap, age increments;
  - acceptance: no boundary object redraw unless topology actually changes.

Exit criteria:

- `tectonics.boundary_objects` persist across archive frames; initial overlap
  matching implemented;
- subduction polarity is recorded; implemented for trench/subduction/active
  margin objects;
- transform/fracture-zone candidates exist as first-class objects; initial
  transform objects carry shear metadata.

Initial implementation:

- `TectonicsModule._boundary_objects` now records stable object ids by matching
  same-kind boundary components against previous boundary objects using cell
  overlap.
- Boundary objects now expose `age_myr`, `persistence`,
  `boundary_continental_fraction`, adjacent `parent_plate_ids`, mean oceanic
  age by plate, continental fraction by plate, and relative-motion metadata.
- Subduction polarity is resolved from local crust context:
  ocean-continent convergence selects the oceanic plate as subducting; ocean-
  ocean convergence selects the older oceanic lithosphere as subducting.
- Transform objects expose `shear_dominance` and are benchmarked with offset
  ridge geometry.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R3 --out out_bench_r3_20260623`
- Latest R3 metrics:
  - `R3.subduction_polarity_old_ocean`: west-old selects plate `0`,
    east-old selects plate `1`;
  - `R3.ocean_continent_margin`: trench continental fraction `0.0`, active
    margin continental fraction `1.0`, polarity basis
    `oceanic_plate_subducts_beneath_continent`;
  - `R3.transform_from_offset_ridge`: transform cell count `12`, lon span
    `41.93 deg`, lat span `8.51 deg`;
  - `R3.boundary_persistence`: stable ridge id persists and age increments
    `25.0 Myr`.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py tests/test_engine.py::test_r3_boundary_objects_expose_persistence_and_polarity_metadata tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> `6 passed in 12.42s`.
- Full suite:
  `.venv/bin/python -m pytest -q`
  -> `45 passed in 397.64s`.

## R4. Ocean Basin And Wilson Cycle Lifecycle

Status: initial ocean-basin lifecycle objects implemented 2026-06-23

Purpose:

- replace current-frame Wilson labels with persistent ocean-basin lifecycles.

Objects:

- `tectonics.ocean_basins`
- `tectonics.rift_systems`
- `tectonics.passive_margins`
- `tectonics.spreading_centers`
- `tectonics.closing_margins`
- `tectonics.sutures`

Lifecycle:

1. rift initiation;
2. narrow sea;
3. spreading ocean;
4. mature ocean;
5. subduction initiation;
6. closing ocean;
7. collision;
8. suture / remnant basin.

Microbenchmarks:

- `R4.rift_to_ocean`
  - fixture: extensional stress across weak continental belt;
  - metrics: rift object births, evolves to narrow sea, then ridge-linked
    ocean basin;
  - acceptance: passive margins attach to both sides and persist.
- `R4.basin_maturation`
  - fixture: stable spreading for several steps;
  - metrics: basin area expands, oceanic age bands grow away from ridge;
  - acceptance: isochron age increases monotonically with distance and time.
- `R4.closure_to_suture`
  - fixture: ocean basin with active trench consuming oceanic crust;
  - metrics: basin area shrinks, gateway restricts, suture forms after
    continent-continent contact;
  - acceptance: basin object does not disappear before closure event.
- `R4.gateway_causality`
  - fixture: narrow seaway between basins;
  - metrics: gateway object has parent basin/margin/rift cause;
  - acceptance: no final terrain-only gateway without parent object.

Exit criteria:

- Wilson-cycle phase field is derived from persistent objects; done via
  `tectonics.ocean_basins` -> `tectonics.wilson_cycles`;
- ocean gateways have parent basin and boundary ids; done;
- current P14 terrain seaway opening is marked as fallback only; still an open
  migration note for R6/R8 terrain integration.

Initial implementation:

- Added lifecycle object sets: `tectonics.ocean_basins`,
  `tectonics.rift_systems`, `tectonics.passive_margins`,
  `tectonics.spreading_centers`, `tectonics.closing_margins`, and
  `tectonics.sutures`.
- Wilson cycles and tectonic ocean gateways are now derived from ocean-basin
  lifecycle objects rather than only from the current-frame boundary kind.
- Lifecycle ids are keyed by plate-pair lineage so opening, closing, and suture
  stages can preserve a basin id across boundary-state transitions.
- Gateways now carry `parent_basin_id`, `parent_boundary_object_id`, and parent
  rift/margin/closing/suture ids.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R4 --out out_bench_r4_20260623`
- Latest R4 metrics:
  - `R4.rift_to_ocean`: basin id persisted as `ocean_basin:plates:0-1`;
  - `R4.basin_maturation`: ridge age p90 `0.0 Myr`, far-ocean age p50
    `125.0 Myr`, basin lifecycle age `25.0 Myr`;
  - `R4.closure_to_suture`: opening, closing, and suture stages keep the same
    basin id;
  - `R4.gateway_causality`: gateway has parent basin
    `ocean_basin:plates:0-1`, parent boundary `passive:gateway`, and parent
    margin.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py tests/test_engine.py::test_r4_wilson_lifecycle_objects_expose_parent_causality tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> `7 passed in 21.17s`.
- Full suite:
  `.venv/bin/python -m pytest -q`
  -> `47 passed in 448.19s`.

## R5. Continental And Margin Object Evolution

Status: initial causal conservation and microbenchmarks implemented
2026-06-23

Purpose:

- replace target-area continental maintenance with object-driven continental
  growth, rifting, accretion, and stabilization.

Objects:

- craton;
- shield;
- platform;
- mobile belt;
- accreted terrane;
- island arc;
- rifted margin;
- suture belt;
- foreland basin;
- large igneous province.

Process model:

- cratonization from old, cool, low-strain continental cores;
- rifting from high rift potential plus inherited weakness;
- accretion from island arcs / terranes reaching active margins;
- collision from closing basin margins;
- old orogens decay into subdued uplands and platforms;
- continental area is diagnostic, not a direct per-step target repair.

Microbenchmarks:

- `R5.craton_survival`
  - fixture: old stable continent under weak boundary noise;
  - metrics: craton area and ids persist; margins change more than cores;
  - acceptance: craton cells are not eroded by area balancing.
- `R5.arc_accretion`
  - fixture: island arc approaches continental margin;
  - metrics: arc terrane attaches, receives parent continent id, changes domain;
  - acceptance: terrane area appears only after collision/accretion event.
- `R5.rifted_margin_split`
  - fixture: supercontinent with high rift potential;
  - metrics: rift opens along weak zone, passive margins form on both sides;
  - acceptance: breakup location follows rift potential, not random seed.
- `R5.no_arbitrary_area_flip`
  - fixture: continental area below target but no accretion/rift/collision
    process active;
  - metrics: no ocean cells become continent solely to hit target;
  - acceptance: land fraction warning is allowed; causality violation is not.

Exit criteria:

- `_conserve_continental` is reduced to a bounded fallback or removed from
  Earth-like primary path; initial implementation now blocks unforced
  ocean-to-continent and continent-to-ocean area repair and records diagnostic
  block flags;
- continent and margin objects explain most coastline changes; partially
  implemented through R1 rift potential, R3/R4 boundary/lifecycle objects,
  craton survival, and active-margin arc accretion.  Expanded shield/platform,
  mobile-belt, foreland-basin, and old-orogen object taxonomy remains for the
  P16/R6 object-derived terrain path.

Initial implementation:

- Added deterministic rifted-margin candidate selection from
  `tectonics.rift_potential` plus inherited lithospheric weakness.
- `_conserve_continental` no longer creates continental crust solely to restore
  target area when no accretion/collision process is active.
- Excess continental area is eroded only near ridge/separating/rifted-margin
  process sources, and craton/stable cores remain protected.
- Balanced shape-aware maintenance is now scoped by process zones when called
  after active conservation.
- Cratonization now requires broad continental interior, sufficient cooling
  age, and no recent reworking, restoring stable cratonic crust in the P12
  Earth-like smoke run without changing the release gate.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R5 --out out_bench_r5_20260623`
- Latest R5 metrics:
  - `R5.craton_survival`: `43` weak-margin cells removed, `0` craton cells
    removed;
  - `R5.arc_accretion`: `43` gained continental cells, `0` outside the
    accretion zone;
  - `R5.rifted_margin_split`: `42`/`42` candidate cells in the imposed weak
    rift belt;
  - `R5.no_arbitrary_area_flip`: target fraction left unmet when no process
    source exists; `0` unforced gained/lost cells.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_r5_continental_margin_evolution_microbenchmarks_pass tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet -q`
  -> passed in targeted runs.

## R6. Terrain From Tectonic Objects

Status: initial object-derived terrain masks and microbenchmarks implemented
2026-06-23

Purpose:

- make elevation, bathymetry, shelves, mountain belts, and hex terrain derive
  from tectonic objects instead of final cleanup rules.

Terrain derivations:

- ridge object -> continuous mid-ocean ridge high;
- transform object -> fracture-zone lineament and ridge offset;
- trench object -> narrow deep trench, forearc, volcanic arc;
- passive margin object -> shelf / slope / rise / sediment wedge;
- collision object -> orogen, plateau, foreland basin;
- rift object -> graben, shoulder uplift, narrow sea;
- craton/platform object -> low-relief interior;
- old orogen object -> subdued upland.

Microbenchmarks:

- `R6.passive_margin_profile`
  - fixture: mature passive margin;
  - metrics: shelf width, slope depth, rise transition, sediment wedge;
  - acceptance: depth increases offshore, no superdeep water at first
    nearshore cell.
- `R6.active_margin_profile`
  - fixture: ocean-continent subduction;
  - metrics: offshore trench, forearc low, volcanic arc high, backarc optional;
  - acceptance: trench is localized and arc sits landward of trench.
- `R6.ridge_transform_bathymetry`
  - fixture: spreading center with transform offsets;
  - metrics: continuous ridge high, fracture-zone lineaments, abyssal plains;
  - acceptance: no dotted ridge province; ridge geometry is object-continuous.
- `R6.orogen_width_and_decay`
  - fixture: collision followed by inactive phase;
  - metrics: orogen width broadens during collision, relief decays after
    activity ends;
  - acceptance: mountain belts are not one-cell checkerboard lines.
- `R6.compiler_consistency`
  - fixture: compiled hex map from all terrain fixtures;
  - metrics: source land/water contradiction, shelf-as-deep-ocean, elevation
    sign mismatch;
  - acceptance: same or stricter than current compiler release gate.

Exit criteria:

- P14/P15 terrain cleanup paths are no longer needed for normal Earth-like run;
  still open.  Initial R6 terrain derivation now consumes persistent
  tectonic objects, but the late Earth-like cleanup passes remain active;
- `elevation.png`, `ocean_depth_provinces.png`, and `hexmap.png` all trace to
  object-derived provinces; partially implemented through object-derived
  ridge/trench/passive-margin/active-margin/transform/suture masks.

Initial implementation:

- Terrain process masks now merge current boundary arrays with persistent
  `tectonics.boundary_objects` and R4/R5 lifecycle objects such as
  `tectonics.spreading_centers`, `tectonics.passive_margins`,
  `tectonics.closing_margins`, `tectonics.rift_systems`, and
  `tectonics.sutures`.
- `_tectonic_relief`, `_regionalize_ocean_floor`, `_ocean_geography`, and the
  open-ocean shoal clamp now consume those process masks instead of reading
  only current-frame boundary arrays.
- Active-margin terrain now produces stronger landward arc relief from active
  margin objects; transform objects create deeper fracture-zone lineaments in
  ocean bathymetry.
- Orogen terrain now distinguishes recent collision belts from old subdued
  orogens: recent orogens get broad uplift, while old orogens keep lower
  residual relief.
- Latest command:
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R6 --out out_bench_r6_20260623`
- Latest R6 metrics:
  - `R6.passive_margin_profile`: shelf p75 `660 m`, slope p50 `1847 m`, rise
    p50 `3164 m`, far-ocean p50 `5196 m`;
  - `R6.active_margin_profile`: trench median depth `3559 m`, landward arc
    relief p75 `324 m`, active-margin classification fraction `0.985`;
  - `R6.ridge_transform_bathymetry`: ridge median depth `2600 m`, transform
    median depth `3600 m`, far abyss median depth `4327 m`;
  - `R6.orogen_width_and_decay`: recent uplift median `608 m`, old uplift
    median `138 m`, broad recent-orogen fraction `0.693`;
  - `R6.compiler_consistency`: source/compiled terrain envelope passed.
- Tests:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_r6_object_derived_terrain_microbenchmarks_pass tests/test_engine.py::test_p12_release_gate_writes_summary_and_contact_sheet tests/test_engine.py::test_earthlike_map_reasonableness_regression -q`
  -> passed in targeted runs.

## R7. Parameter Calibration Matrix

Purpose:

- tune model coefficients through microbenchmarks and small release matrices.

Parameter groups:

- mantle heat diffusion and cooling;
- continental insulation coefficient;
- slab-pull coefficient;
- ridge-push coefficient;
- basal-drag coefficient;
- collision resistance;
- subduction initiation threshold;
- rift initiation threshold;
- cratonization age/stability thresholds;
- passive-margin sedimentation rate;
- transform segmentation scale;
- erosion and relief-decay rates.

Microbenchmarks:

- `R7.parameter_sensitivity`
  - fixture: run selected R1-R6 benchmarks over parameter grid;
  - metrics: monotonicity, stability, and geological plausibility bands;
  - acceptance: no parameter has hidden inverted response in its primary
    benchmark.
- `R7.multi_seed_ensemble`
  - fixture: Earth-like short ensemble with deterministic process rules and
    small allowed perturbations;
  - metrics: land fraction, component count, largest component, ridge/trench
    ratio, passive/active margin ratio, ocean age distribution;
  - acceptance: variance is explainable by perturbations, not random event
    placement.
- `R7.baseline_earth_analogue`
  - fixture: 8000-cell Earth-like run;
  - metrics: release-gate fields plus refactor-specific object causality
    metrics;
  - acceptance: no hard regressions versus current P15 baseline.

Exit criteria:

- tuning values are stored in preset/spec parameters or an explicit parameter
  object;
- benchmark summaries include parameter hashes for reproducibility.

## R8. Release Integration And Migration

Purpose:

- move the refactored process path from experimental to default without losing
  existing stable outputs.

Migration strategy:

- introduce experimental flags first, for example `tectonics_refactor_level`;
- keep current P15 path as fallback until R1-R6 benchmarks pass;
- replace one major-event source at a time;
- compare old and new paths with release-gate summaries and benchmark reports.

Microbenchmarks:

- `R8.old_new_comparison`
  - fixture: current Earth-like seed 42 plus two alternate seeds;
  - metrics: release status, object causality, morphology, ocean geography,
    archive continuity;
  - acceptance: refactored path matches or improves hard-gate behavior.
- `R8.archive_causality`
  - fixture: multi-frame archive from R4/R5 lifecycle worlds;
  - metrics: percentage of terrain provinces with parent object/event;
  - acceptance: parentless major terrain provinces below threshold.
- `R8.earth_geomorphology_coverage`
  - fixture: run the E1-E5 Earth coverage suites and one Earth-like reference
    world;
  - metrics: covered, partial, weak, and missing feature counts; parentless
    major landform fraction; ridge-transform continuity; passive/active margin
    profile scores; ocean-age isochron score; orogen lifecycle score;
  - acceptance: no non-ice required feature remains `none`, all plate-
    boundary and ocean-basin features reach at least `partial`, and generated
    parentless major landform area remains below the current hard gate;
  - status: initial implementation writes `earth_geomorphology_coverage` to
    P12 release summaries; generated parent-object linkage has reduced the
    2500-cell Earth-like `parentless_major_landform_fraction` to `0.037`.
- `R8.performance_budget`
  - fixture: benchmark and 8000-cell Earth-like run;
  - metrics: runtime and memory;
  - acceptance: microbenchmark suite remains fast enough for routine tuning,
    and full release remains practical.

Exit criteria:

- default Earth-like run uses principle-model major-event decisions;
- current map assets are generated from persistent tectonic objects;
- `earth_geomorphology_coverage` is written to the release summary;
- release gate and microbenchmark summaries are both part of development
  review.

## Cross-Stage Benchmark Summary

| Stage | Main Replacement | Key Microbenchmark Signal |
| --- | --- | --- |
| R0 | hidden randomness -> audited randomness | zero unclassified RNG calls |
| R1 | random plume/rift locations -> potential fields | plume/rift at field maxima |
| R2 | random pole refresh -> torque proxy | velocity responds to slab/ridge/collision terms |
| R3 | cell boundary labels -> persistent boundary objects | stable ids, polarity, transform continuity |
| R4 | frame Wilson labels -> basin lifecycle | rift-to-ocean-to-suture object continuity |
| R5 | area flip conservation -> object evolution | no continent creation without process cause |
| R6 | final terrain cleanup -> object-derived terrain | ridge/trench/margin/orogen profiles pass |
| R7 | ad hoc tuning -> parameter calibration | monotonic response and plausibility ranges |
| R8 | experimental path -> default path | release gate, causality, and Earth coverage benchmarks pass |

## Earth Coverage Benchmark Summary

| Suite | Feature focus | Key signal |
| --- | --- | --- |
| E1 | craton, platform, interior basin, old orogen, foreland basin | implemented; continental interiors have object causes and lifecycle states |
| E2 | passive margin, active margin, shelf, delta/fan | implemented; margin profiles match tectonic setting and sediment supply |
| E3 | ridge, transform, fracture zone, abyssal plain, ocean age | implemented; ocean-basin fabric is continuous and history-derived |
| E4 | island arc, backarc, terrane, plume, hotspot track, LIP | implemented; arc/plume features follow process fields and plate motion |
| E5 | ice sheet, glacial erosion, rebound | implemented; cryosphere landforms appear only where ice/surface-process fields allow |

## Immediate Next Work

R0 through the initial R6 object-derived terrain path are implemented.  E1
continental-interior, E2 margin/shelf, E3 ocean-basin fabric, E4 arc/plume,
E5 cryosphere/surface-process geomorphology coverage, and the initial R8
`earth_geomorphology_coverage` release-summary integration with tightened
generated parent-object linkage are also implemented and passing.  R8
generated ocean-fabric fallback now makes `transform_fault` and
`fracture_zone` visible in Earth-like worlds that do not yet emit explicit
transform boundary objects.  R8 continental fallback now makes
`interior_basin`, `foreland_basin`, and `old_subdued_orogen` visible in
generated Earth-like worlds from stable platform depocenters, recent orogenic
loading, and eroded old-orogen envelopes.
P16/R8 now makes `plateau` visible in generated Earth-like worlds from broad
thickened collision/LIP highland cores.  P16 now also has deterministic
continental shape-pressure maintenance and a seam-safe connected-ocean
major-basin partition.  P16/P18 final-step scheduler refresh now keeps the
final world time synchronized with final tectonics, terrain, and object state,
and process-belt continental detail classification removes the prior
orogen/plateau overpaint warning.  P19 fixed background shape-maintenance
margin fills so they inherit parent continental properties instead of being
mislabelled as new arcs; active accretion remains arc/suture-like.  The P19
audit also ruled out two tempting but wrong fixes for the remaining ribbon
warning: broad drowning of narrow land and aggressive shallow-water platform
expansion both worsen morphology in counterfactual checks.  The 2500-cell
Earth-like P12 gate is `warn`
with no hard failures, parentless major landform area at `0.0`, multiple major
ocean basins, and `orogen_or_plateau_fraction_of_land=0.354`.  P20 now exposes
resolved plate topology objects and uses them to guide deterministic
large-plate breakup from topology pressure, rift potential, boundary load, and
crustal weakness.  This removed the stochastic local plate reorganization split
seed/pole/rate path and lowered the 2500-cell Earth-like land ribbon fraction
from the P19 baseline `0.562` to `0.399`, though ribbon warnings remain.  P20
also now records topology-aware microplate capture metadata and distinguishes
oceanic microplate capture from continental-cargo capture.  The important
engineering finding is that tiny microcontinents should not reserve scarce
plate ids at this stage: doing so blocks large-plate breakup and worsens
morphology.  The current policy captures the plate label while preserving the
crust/continent cargo for the continent/terrane object layer.  P20 now also
adds deterministic continent/terrane lifecycle events from id-overlap and
capture metadata: `continent_split`, `continent_merge`, `continent_birth`,
`continent_loss`, `terrane_capture`, and `microcontinent_plate_capture`.
The latest P20 pass also adds deterministic passive-margin/platform
progradation and repairs generated-world `interior_basin` object ownership so
detail-level basin cells are not preempted by old-orogen objects.  This keeps
continental conservation and non-ice landform visibility stable at the 2500
cell regression scale.  P21 now adds `tectonics.breakup_seaways`, a
deterministic supercontinent-breakup / multi-rift seaway object layer, and a
terrain response that only opens object-backed corridors when they reduce
largest-landmass dominance.  The second P21 pass extends accepted breakup
objects to coast/boundary paths through weak cells, regionalizes continental
sediment into rift/foreland/margin depocenters, and relaxes overbroad mature
suture/accretionary collages so old collision belts no longer turn most land
into high plateau.  The latest 8000-cell Earth-like routine review is `warn`
with `0` failed entries: land fraction `0.235`, component count `13`, ribbon
fraction `0.549`, and ocean basins `6`.  Ribbon and coastline complexity
remain warning-level, so P21 is not complete as a morphology phase.  The first
24000-cell medium audit has now been run and failed: land fraction `0.238`,
component count `11`, largest land component `0.944`, ribbon fraction `0.482`,
and ocean basins `26`.  The audit proves that some ribbon artifacts improve
with resolution, but the main landmass remains a single connected
supercontinent.  A direct terrain-layer completion attempt did not change the
failure and should not be retained as the solution.  P22 now has a passing
topology-scored benchmark for this pattern: candidate breakup objects record
expected partition topology, generate interior rift-axis candidates, and reject
oversized weak blankets.  The P22 fixture now selects the interior rift corridor
(`topology_score=1.029`, object split balance `0.811`), terrain opens one
object-backed seaway, the largest exposed-land component falls to `0.537`, and
interior-rift opening reaches `1.000`.  A fresh P22 24000-cell audit has also
been run.  It improves the largest component from the P21 audit's `0.944` to
`0.856`, but still fails the release gate as an effective single
supercontinent; ribbon remains high at `0.487`, mean land elevation is
`2131 m`, and p95 land elevation is `5585 m`.  The next phase is therefore P23:
generated-world residual supercontinent repair, driven by explicit
breakup/rift/seaway object telemetry rather than subjective map inspection.
P23a/P23b now add focused coverage for medium rifted components, continental
divergent-boundary lifecycle, and boundary-seeded medium breakup objects.  On
the real 24000 generated world this raises rift-system objects from `0` to
`12`, but they remain skeletal (`3-4` cells each); breakup-seaway count is
still `1`, seaway area remains about `0.0037`, and largest exposed-land
component is still `0.852`.  P23c now adds component-level telemetry, a
multi-corridor breakup fixture, and a guarded rift-corridor apron expansion
that keeps P22/P23 microbenchmarks green.  The generated 24000 telemetry,
however, still has only one accepted breakup seaway: the dominant component has
`10` candidates but only `1` accepted object, terrain opens only `0.0036` of
global area, and the largest exposed-land component remains `0.852`.  P24
should therefore target the object-to-terrain conversion and exposed-land
connectivity check, not higher resolution alone.  P24a now fixes and tests one
real conversion failure mode: a breakup/rift axis that has already become ocean
must still act as a process source for seaway propagation into adjacent weak
land.  The P24 fixture passes, but the 24000 generated world reports
`terrain_breakup_seaway_source_reuse=0` and largest exposed-land component
remains `0.852`, so the current seed is not blocked by submerged-axis reuse.
P24b adds per-object terrain attempt telemetry and stage-level land-component
telemetry; that evidence showed the accepted breakup corridor worked
(`largest_share 0.945 -> 0.476`) but `_regionalize_ocean_floor` later raised
the opened shallow seaway back into a land bridge.  P24b now records opened
object-backed corridors and protects the corridor plus one ocean-neighbour
apron during ocean-floor regionalization.  P25a then calibrates continental
hypsometry and highland overpaint: mature suture/LIP/accretionary collages and
stable cratons are subdued unless active process evidence supports high relief,
and suture/LIP land no longer becomes highland province by default.  The latest
24000 P12 audit is `warn` with `0` failed entries: land fraction `0.237`,
component count `13`, largest land component `0.483`, ribbon `0.490`, mean land
elevation `1248 m`, p95 `3116 m`, high land above `2500 m` `9.8%`, and
suture/LIP/highland province share `9.3%`.  The residual work is no longer a
single-supercontinent or false-highland hard failure; P25b should target
land-fraction undershoot, ribbon/coastline complexity, a small mean-elevation
excess, and slightly high abyss fraction.
Non-ice fixture-only feature classes are now generated-world visible in the
2500-cell Earth-like gate; the remaining fixture-only classes are the deferred E5
cryosphere landforms.  The next concrete task should reduce the remaining
upstream morphology warnings surfaced by that gate:

1. keep `E1.craton_platform_basin`, `E1.old_orogen_decay`,
   `E1.collision_plateau`, and `E1.foreland_basin` in regression; current suite
   `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E1 --out out_bench_e1_20260623`
   passes `4/4`;
2. keep `E2.passive_margin_shelf_wedge`, `E2.active_margin_trench_arc`, and
   `E2.delta_fan` in regression; current suite
   `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E2 --out out_bench_e2_20260623`
   passes `3/3`;
3. keep `E3.ridge_transform_fracture_zone`,
   `E3.abyssal_plain_sedimentation`, and `E3.ocean_age_isochrons` in
   regression; current suite
   `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E3 --out out_bench_e3_20260623`
   passes `3/3`;
4. keep `E4.island_arc_accretion`, `E4.back_arc_basin`,
   `E4.hotspot_track`, and `E4.large_igneous_province` in regression; current
   suite
   `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E4 --out out_bench_e4_20260623`
   passes `4/4`;
5. keep `E5.ice_sheet_loading`, `E5.glacial_erosion`, and
   `E5.postglacial_rebound` in regression; current suite
   `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E5 --out out_bench_e5_20260623`
   passes `3/3`;
6. keep generated-world visibility for non-ice feature classes in regression;
   keep E5 cryosphere generated-world expression deferred until climate/ice
   work resumes;
7. keep generated `parentless_major_landform_fraction` below `0.10` while
   tightening toward zero;
8. keep `P21.breakup_seaway_objects` and
   `P21.object_driven_terrain_seaway` in regression; current suite passes and
   proves the breakup corridor is object-backed rather than a renderer cut;
9. keep the expanded P21 suite in regression: breakup object generation,
   object-driven terrain seaway response, and continental sediment depocenter
   contrast now pass; the 24000-cell P21 medium morphology audit has been run
   and failed on a single-supercontinent hard gate, which defined P22;
10. keep the new P22 topology-scored supercontinent breakup regression green:
   interior rift connectivity must outrank peripheral weak-neck cuts when it is
   the only candidate that partitions the main continent;
11. keep P12 `tectonic_object_telemetry` in the release summary so 8000/24000/
   72000 audits can distinguish missing breakup objects from terrain non-
   opening; the field must include rift-system count, breakup-seaway count,
   topology/split metrics, terrain seaway openings, and area fractions;
12. keep P24 focused on the residual 24000 failure after P23c: breakup objects
   can be topology-capable in continental-crust space but still fail to sever
   exposed-land connectivity because terrain opening is too narrow or because
   adjacent non-continental/terrane bridges keep the landmass connected; P24a
   covers submerged-axis source reuse, while P24b exposes terrain-candidate
   effectiveness and fixes ocean-floor regionalization reconnecting opened
   object-backed seaways;
13. continue migrating remaining P14/P15 coastline and ribbon cleanup from late
   terrain patches into continent, margin, rift, platform, craton, old-orogen,
   and ocean-basin object responses; the current hard supercontinent,
   single-basin, final-state-sync, detail-overpaint, stochastic large-plate
   split, unlabeled microplate-capture, and missing continent/terrane lineage
   problems are fixed, but ribbon warnings remain; the next implementation step
   should validate the P20 object model at 8000/24000 cells and then improve
   continent/margin object geometry, not add another final terrain cleanup pass;
14. use the implemented `profile-resolution` command as the resolution-ladder
   preflight for geomorphology review: keep 900-2500 cell
   microbenchmarks for mechanism tests, 8000 cell Earth-like runs for routine
   release review, 24000 cell medium morphology audits after major object
   changes, and occasional 72000 cell offline deployment tests for small-scale
   features such as isthmuses, deltas/fans, short mountain systems, narrow
   shelves, straits, and island-arc spacing after their parent process objects
   exist at lower resolution; the summary now records the 8000/24000/72000
   projection targets and the optional acceleration packages visible in the
   environment; 72000-cell runs are useful for small-feature visibility and
   performance scaling, but they must not be used to claim macro-topology
   quality while the 24000 largest-landmass gate still fails; each 72000 audit
   should compare the same seed at 8000/24000/72000, start with core generation
   plus diagnostics before rendering, and attach parent-object provenance for
   every claimed small-scale feature;
15. profile and optimize the CPU path before adding optional acceleration; the
   new `resolution_profile_summary.json` output should identify whether
   repeated graph width/component passes, terrain smoothing, hydrology routing,
   `cKDTree` plate rasterization, map compilation, or rendering is the first
   scaling target for the chosen cell count; Numba should be considered first
   for graph-heavy loops or edge-reduction hot paths, while CuPy/JAX should
   only be added behind explicit backend boundaries for dense field kernels
   after numerical parity and end-to-end runtime benefit are measured;
16. keep climate, ocean-current, and monsoon redesign paused except for
   regression fixes needed to keep existing tests passing.

P20 should use professional-tool and real-data anchors:
GPlates/pyGPlates-style finite rotations, plate circuits, topological features,
and resolved topology objects for plate/continent representation
(`https://www.gplates.org/docs/pygplates/`); ETOPO-derived distributions for
hypsometry, shelves, abyssal plains, mountain belts, and coastline complexity
(`https://www.ncei.noaa.gov/products/etopo-global-relief-model`,
DOI `10.25921/fd45-gt74`).  The goal is not exact Earth replay, but
Earth-plausible feature distributions with persistent causal objects and
parameter-ledger tuning.  Initial P20 implementation now writes
`tectonics.plate_topologies`, a resolved topology summary per active plate with
area, component count, adjacent plate ids, shared-edge counts, boundary-load
fractions, continental/terrane cargo, motion source, and R2 force components.
The second P20 step now consumes these objects for deterministic large-plate
split decisions, using topology pressure, rift potential, boundary load, crust
type, and stability instead of random split seed/pole/rate jitter.  The third
P20 step adds topology-aware microplate capture metadata and a cargo policy:
tiny plate labels can be captured by dominant neighbouring plates while
continental or terrane cargo remains in crust/object ids for later lifecycle
storytelling.  The fourth P20 step adds those lifecycle records from
deterministic id overlap and capture metadata.  The first `8000`-cell
Earth-reference validation has now been run and it failed on land ribbons:
exposed-land ribbon fraction `0.677`, narrow necks `22.52` per 1000 land
cells, largest-coastline complexity `16.48`, exposed-land width p50/p90 `2/4`
steps, and continental-crust width p50/p90 `3/5` steps.  The next P20 steps
are therefore targeted continent/margin geometry changes, direct ETOPO raster
sampling, then 24000 validation before any optional 72000-cell offline audit.
Initial ETOPO-anchored screening is now present in P12 as
`earth_reference_distribution`; at 2500 cells it shows the strongest current
Earth-reference deviations are continental: low land fraction, high land
ribbon fraction, excessive mean land elevation, too much land above 2500 m, and
slightly high orogen/plateau coverage.  At 8000 cells the continental geometry
problem worsens while ocean province distribution remains inside the initial
broad envelope, so the next implementation target is upstream continent,
margin, platform, rift, and accretion object geometry rather than a final
terrain cleanup pass.

P25a and P25b have now narrowed the generated-world mismatch.  P25a suppresses
broad false mature highlands and brings 24000-cell land hypsometry close to the
screening envelope.  P25b adds seaway-safe coastline land-budget payback:
modern coastline smoothing can remove unsupported island-chain area, then pay
that area back into same-component shallow continental shelves without
refilling object-backed breakup seaways.  The latest 24000-cell Earth-like
audit is `warn`, not release-clean: land fraction is now inside the screening
range at `0.253`, exposed land components are down to `4`, largest component
is `0.486`, and mean land elevation is `1134 m`; remaining
Earth-reference misses are exposed-land ribbon fraction `0.410`,
largest-landmass coastline complexity `21.46`, and abyss fraction `0.715`.
The next refactor step should be P26: widen and simplify continent/margin
geometry through process objects and parameter calibration, and tune ocean
province allocation, rather than adding another area-only terrain cleanup.
An initial P26 spike that widened passive-margin progradation/shape-maintenance
passed isolated fixtures but failed integration, so it was reverted.  The
failed 24000 audit (`out_p26_earthlike_24000_upstream_margin_geometry_20260624`)
had land `0.216`, components `16`, ribbon `0.565`, coastline complexity
`34.53`, and `40` basins.  The restore audit
(`out_p25b_restore_earthlike_24000_after_p26_revert_20260624`) returns to the
accepted P25b baseline: land `0.253`, components `4`, ribbon `0.410`, coastline
complexity `21.46`, basins `14`, and `0` failures.  Therefore P26 must be
designed with an integration gate from the start: no isolated microbenchmark
can count as accepted unless the 8000/24000 Earth-like audits do not regress
P25b land fraction, component count, ribbon, basin count, or seam continuity.
The integration gate is now codified as
`aevum.diagnostics.p26_regression_gate`; it compares P12 summaries and rejects
candidate runs that regress release/validation status, land fraction,
component count, land/continental ribbon, largest-coastline complexity,
ocean-basin count, or basin seam continuity.  A P26 candidate should be
documented only after the 24000-cell P12 summary passes this comparator
against
`out_p25b_restore_earthlike_24000_after_p26_revert_20260624/p12_tectonics_release_summary.json`.
P12 summaries now also expose `p26_ribbon_drivers`, a read-only component
attribution block that ranks the exposed-land and continental-crust components
responsible for ribbon area and reports their domain/origin/detail/stability
mix.  The next P26 implementation should begin from that attribution rather
than applying another uniform margin-widening or final cleanup pass.
The first 8000-cell attribution run
(`out_p26_driver_audit_earthlike_8000_20260624`) identifies young
suture/accreted-terrane continental crust as the dominant residual driver.
Two local mechanism probes were rejected after world-level audits: broadened
quiet-collage maturation worsened ribbon to `0.465`, and simple boundary
rework localization failed hard with ribbon `0.575`.  Those behavior changes
were removed.  `p26_ribbon_drivers` now includes time-since-rework,
recent-orogeny, recent-volcanism, active-rework, and quiet inherited
arc/suture shares, plus summary-level temporal driver hints.  The rerun
`out_p26_time_attribution_earthlike_8000_20260624` shows active recent
collision/subduction rework, not quiet inherited collage, is the dominant
current driver: continental ribbon active-rework share is about `0.91` while
quiet inherited arc/suture share is about `0.00`.  The next P26 step should
therefore design a microbenchmark for over-broad recent rework swaths before
attempting another production mechanism change.
That microbenchmark now exists as the diagnostic-only `P26` suite.  It uses
`aevum.diagnostics.p26_rework_footprint` to flag broad recent rework outside a
boundary corridor and to accept localized active belts.  P12 summaries also
include `p26_rework_footprint`, so the next production change can be evaluated
both in fixture space and in Earth-like world-level audits.
Two production probes using a boundary-corridor rework core were then rejected:
one improved exposed-land ribbon but regressed continental-crust fragmentation,
continental ribbon, exposed oceanic land, and craton fractions; the other
failed the 8000-cell audit outright.  The next production design should
therefore introduce an explicit deforming-network or province state that can
separate active provenance from continental topology, instead of only narrowing
the origin/rework write mask.
That state layer has now been added as `tectonics.deformation_intensity`,
`tectonics.deformation_style`, and `tectonics.deforming_networks`.  It is not
yet used to alter terrain or crust provenance; it exists so the next P26
production change can read active deformation directly without fragmenting the
continental-crust topology through `crust.origin`.  The P12 summary now exposes
`p26_deforming_networks` telemetry so candidate production changes can be
screened for broad active deformation overprint, core/shoulder balance, and
ribbon coupling before a 24000-cell audit is treated as meaningful.  The
deformation axes are now thinned from plate-contact networks so terrain does
not amplify raw rasterized separating swaths into broad rift relief.  The first
terrain-side consumer is limited to province/detail/object semantics; direct
elevation response is still pending because it must not lower land fraction,
increase exposed-land components, or over-fragment ocean basins in the P26
regression gate.  The first accepted elevation-side step is therefore a small
constrained relief correction only on already exposed broad continental
deformation cores/shoulders, with stage telemetry proving no land-mask change
and the same-configuration 24000-cell P26 gate passing.  The accepted P26
version now also requires interior exposed-land width and pre-existing
elevation support, so active but narrow island/bridge cells and active but
low-elevation continental cells are negative controls in the P26 microbenchmark
suite.
The next accepted step, P27, adds explicit terrain response state:
`terrain.orogenic_load` and `terrain.foreland_accommodation`.  This is the
first object-level terrain-province bridge after P26: continental detail and
`foreland_basin` objects consume one orogen-to-foreland response state and keep
parent tectonic ids.  Direct sediment/elevation feedback from that state was
tested and rejected at 24000 cells because it regressed Earthlike land fraction;
future sediment coupling needs a separate conservation gate before it can be
accepted.
P28 now provides the sediment-coupling budget gate that must pass before
`foreland_accommodation` can affect production sediment/elevation again.  The
gate is diagnostic-only and currently passes foreland and passive-margin
fixtures with conserved sediment volume and zero projected land-mask change;
production terrain still does not consume it.  P29 should address the visible
inland-flatness gap directly.  The target is not random interior noise: use
ETOPO/GEBCO-style hypsometry references and mechanism fixtures for cratonic
shields, platforms, interior sedimentary basins, old-subdued orogens, rift
valleys, plateau margins, LIP surfaces, and erosional escarpments, while
preserving P26/P27 land-fraction and basin-count regression gates.
The first P29 implementation is now in production behind benchmark gates:
`P29.inland_mechanism_diversity_fixture`, `P29.inland_relief_response_preserves_land_mask`,
`P29.flat_interior_detector`, and `P29.unparented_speckle_detector` distinguish
a mechanism-rich broad continent from a one-elevation flat interior and from
high-relief texture with no tectonic parentage.  Production terrain now consumes
that signal through a constrained inland-relief pass that only edits already
exposed broad continental interiors and clamps depressions above sea level.  The
latest 8000-cell follow-up remains `warn` with `0` hard failures and no accepted
planform regression, but it also confirms that the broad interior platform mode
is still over-compressed: inland p25/p50/p75 elevations cluster near
`1087/1090/1123 m`, with only `36 m` IQR.  A wider old-orogen postprocess was
tested and rejected because it did not improve inland IQR and worsened
ribbon/basin metrics.  The remaining refactor work is therefore upstream:
reduce ribbon land and over-complex coastlines, add persistent continental
interior province/state for shields, sag basins, platform swells, old-orogen
roots, rift shoulders, escarpments, and dynamic topography, improve
generated-world plateau expression, and move P28 sediment coupling from
diagnostic to production only after its mass-balance constraints are preserved.
P30 now keeps that inland-state layer tested, including a moderate
platform-swell fixture distinct from plateau margins.  P31 adds the next
planform-specific gate: `P31.broad_multicontinent_planform_reference`,
`P31.ribbon_land_planform_detector`, and
`P31.overcomplex_coastline_planform_detector`.  The current 24000-cell P29
baseline is inside the initial envelope for land fraction, major landmass
partition, and component count, but remains outside for ribbon fraction
(`0.424` vs max `0.35`) and coastline complexity (`23.276` vs max `8.0`).
P32 has added the first coastline-complexity production fixture plus
object-constrained local simplification guards.  The fixture reduces coastline
complexity `12.63 -> 11.71` and ribbon fraction `0.363 -> 0.329` while
preserving protected seaways, and the 8000-cell generated-world smoke remains
release-gate `warn` with no hard failures.  However the real generated world
still has component count `3`, largest land component `0.648`, ribbon `0.486`,
and largest coastline complexity `15.443`; the local P32 swap telemetry is
still `0.0` in that run.  Therefore the next refactor should move upstream into
P33 ribbon/partition work: reduce ribbon land and over-dominant landmasses
through continent/margin lifecycle mechanisms while preserving seaways, ocean
basins, island arcs with tectonic parents, and antimeridian seam continuity.

Do not rewrite the full plate solver in one pass.  Replace one major random
decision at a time after its benchmark proves the process model in isolation.
