# Earth Geomorphology Coverage Contract

Status: active coverage contract
Owner: tectonics / terrain / diagnostics integration
Created: 2026-06-22

This document defines the Earth-derived geomorphology that Aevum's tectonics
and terrain model must be able to generate.  It is not a requirement to copy
modern Earth.  It is a requirement that the process model can generate the same
classes of landforms for the right geological reasons.

The staged research, source extraction, implementation, and microbenchmark plan
for turning this coverage contract into testable reference metrics is archived
in `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`.

Reference visual material currently archived locally:

- `/Users/rayw/Projects/aevum/out_reference_earth_20260622/nasa_topography_5400.jpg`
- `/Users/rayw/Projects/aevum/out_reference_earth_20260622/nasa_bathymetry_5400.jpg`
- `/Users/rayw/Projects/aevum/out_reference_earth_20260622/earth_real_vs_aevum_p15_comparison.png`

## Coverage Rule

A landform is considered covered only when all of the following are true:

- it can be produced in at least one deterministic microbenchmark fixture;
- it has a parent process object or persistent lifecycle object;
- its diagnostic metrics match a plausible Earth-derived range or relationship;
- its rendered map expression is recognizable at both raster and compiled-hex
  scales;
- it can be explained through archive provenance, not by final-map repair.

## Required Feature Matrix

| Earth feature class | Process model required | Required objects / fields | Current status |
| --- | --- | --- | --- |
| Stable cratons and shields | old continental cores, low strain, cool lithosphere, cratonization | `tectonics.cratons`, crust age, crust stability, terrain shield province | E1 initial |
| Continental platforms | old craton-adjacent interiors, sediment cover, low relief | continent object, platform detail, erosion/sediment state | E1 initial |
| Interior sedimentary basins | flexure, extension, subsidence, sediment routing | basin object, sediment thickness, lowland terrain detail | E1 initial |
| Rift valleys and grabens | extensional stress plus inherited weakness and heat anomaly | rift system object, rift potential, normal-fault terrain | weak |
| Passive margins | rifted margin lifecycle, thermal subsidence, sediment wedge | passive-margin object, shelf/slope/rise fields | E2 initial |
| Active continental margins | ocean-continent subduction, trench, forearc, volcanic arc | subduction object, polarity, trench, arc, active margin | E2 initial |
| Collision orogens | continent-continent convergence, crustal thickening, suture | collision object, suture belt, orogen object | world |
| High plateaus | thickened crust, low erosion or post-collision support | plateau object/province, crust thickness, uplift history | world |
| Old subdued orogens | inactive collision belts, erosion, thermal relaxation | inherited orogen object, relief-decay state | E1 initial |
| Foreland basins | flexural loading beside orogens | foreland basin object, sediment routing | E1 initial |
| Island arcs | oceanic subduction, arc volcanism, possible accretion | island-arc object, volcanic arc, terrane id | E4 initial |
| Back-arc basins | retreating subduction, extension behind arc | back-arc basin object, young oceanic crust | E4 initial |
| Accreted terranes | arc/microcontinent collision and attachment | terrane object with parent continent and event | E4 initial |
| Large igneous provinces | mantle plume or decompression melting | plume potential, plume object, LIP object | E4 initial |
| Hotspot island/seamount chains | plume track plus plate motion | plume track object, age-progressive chain | E4 initial |
| Mid-ocean ridges | divergent boundary, spreading center lifecycle | spreading-center object, ridge terrain, age isochrons | E3 initial |
| Transform faults | ridge offsets and shear-dominant plate motion | transform object, fracture-zone terrain | E3 initial |
| Oceanic fracture zones | inactive transform traces across oceanic crust | fracture-zone object, bathymetric lineament | E3 initial |
| Abyssal plains | old oceanic crust plus sediment smoothing | ocean basin object, abyss province, sediment cover | E3 initial |
| Ocean trenches | subduction polarity and slab pull | trench object, depth province, active-margin link | E2 initial |
| Forearc / accretionary prisms | sediment scraped at subduction margin | forearc/prism object, active margin terrain | E2 initial |
| Continental shelves | passive/active margin subsidence and sedimentation | shelf object/field, shelf width, sediment wedge | E2 initial |
| Submarine fans and deltas | river sediment routing to margin | delta/fan object, sediment discharge | E2 initial |
| Marginal seas and gateways | basin lifecycle, rifting, arcs, changing sea connections | basin/gateway object with parent lifecycle | partial |
| Glacial ice sheets and glacial erosion | climate-ice coupling and bedrock modification | ice object, glacial erosion/deposition fields | E5 initial |

## Current Coverage Summary

The P15 baseline is suitable as a stable geography fixture, but not yet as a
complete Earth-like geomorphology generator.

Covered enough for downstream use:

- broad land/ocean separation;
- cratons / continental interiors / margins at coarse scale;
- initial object-level shields, platforms, interior basins, subdued old
  orogens, and foreland basins in deterministic fixtures;
- initial object-level passive margin sediment wedges, active margin trench /
  forearc / volcanic-arc complexes, and delta/fan objects in deterministic
  fixtures;
- initial object-level spreading centers, transform faults, fracture zones,
  abyssal plains, and ocean-age isochron bands in deterministic fixtures;
- initial object-level island arcs, accreted terranes, back-arc basins,
  age-progressive hotspot tracks, and LIP objects in deterministic fixtures;
- initial object-level ice-sheet loading, glacial erosion, and postglacial
  rebound diagnostics in deterministic fixtures;
- collision and arc signals;
- shelves, slopes, abyssal plains, trenches, ridges at coarse scale;
- terrain/hex compiler consistency.

Weak or incomplete:

- rendered world-scale expression of E1 interior landform objects;
- rendered world-scale expression of E2 margin landform objects;
- rendered world-scale expression of E3 ocean-basin fabric objects;
- rendered world-scale expression of E4 arc/plume landform objects;
- rendered world-scale expression of E5 cryosphere/surface-process landform
  objects;
- persistent ocean-basin lifecycle validation across long generated archives;
- integration of cryosphere-landform objects with the future climate/ice
  redesign.

Temporary fallback that must not count as coverage:

- final terrain seaway cutting;
- final coastline smoothing;
- arbitrary component cap / land payback;
- final open-ocean shoal clamp;
- current-frame Wilson-cycle labels with no persistent basin lifecycle.

## Current Development Priority

Status updated 2026-06-23 after initial E5 cryosphere/surface-process
microbenchmarks.

The next work should prioritize real-Earth geomorphology coverage, not climate,
ocean-current, or monsoon design.  Climate-facing fields should receive only
regression fixes needed to keep existing tests passing until the geomorphology
stack has stronger object coverage.

Use the archived real-Earth topography/bathymetry references and the feature
matrix above to drive implementation order:

1. **Continental interiors**: stable cratons/shields, platforms, interior
   sedimentary basins, old subdued orogens, and foreland basins.  These are
   needed to make large continents read like real Earth rather than smooth
   plates plus boundary mountains.
2. **Continental margins**: passive-margin shelf/slope/rise/sediment wedges,
   active-margin trench/forearc/arc/accretionary-prism geometry, and deltas or
   submarine fans where river sediment reaches margins.
3. **Ocean-basin fabric**: continuous ridge-transform-fracture-zone networks,
   age-derived abyssal plains, trenches tied to polarity, and hotspot/seamount
   tracks tied to plume and plate motion.

The immediate benchmark path is now E1-E5 integrated Earth-like release-gate
coverage.  E1 through E5 have initial implementations and should remain in
regression:

- `E1.craton_platform_basin` - implemented, passing
- `E1.old_orogen_decay` - implemented, passing
- `E1.collision_plateau` - implemented, passing
- `E1.foreland_basin` - implemented, passing
- `E2.passive_margin_shelf_wedge` - implemented, passing
- `E2.active_margin_trench_arc` - implemented, passing
- `E2.delta_fan` - implemented, passing
- `E3.ridge_transform_fracture_zone` - implemented, passing
- `E3.abyssal_plain_sedimentation` - implemented, passing
- `E3.ocean_age_isochrons` - implemented, passing
- `E4.island_arc_accretion` - implemented, passing
- `E4.back_arc_basin` - implemented, passing
- `E4.hotspot_track` - implemented, passing
- `E4.large_igneous_province` - implemented, passing
- `E5.ice_sheet_loading` - implemented, passing
- `E5.glacial_erosion` - implemented, passing
- `E5.postglacial_rebound` - implemented, passing

Climate-system work such as seasonal monsoon derivation, ITCZ migration, storm
tracks, and coupled ocean circulation remains paused except for regression
fixes.

Resolution note:

- The normal 2500-8000 cell Earth-like runs are not sufficient to visually prove
  every small-scale real-Earth landform.  They remain valid for process
  microbenchmarks, major continent/ocean topology, and release-gate regression.
- Features such as narrow isthmuses, small deltas/fans, short mountain systems,
  straits, compact marginal seas, and island-arc spacing should receive
  occasional `24000` and, once performance allows, `72000` cell deployment
  reviews.
- Before a `24000` or `72000` deployment review, run
  `python -m aevum.cli profile-resolution` for the chosen preset/cell ladder
  and archive `resolution_profile_summary.json` with the maps.  This keeps
  visual coverage claims tied to measured build/run/compile/diagnostic costs,
  high-resolution runtime projections, and the optional acceleration packages
  visible in the environment.
- A `72000` review should compare the same seed at `8000`, `24000`, and
  `72000` cells.  It should begin with generation plus diagnostics and only add
  compilation/rendering after the core model timing is acceptable.  Claims
  about isthmuses, deltas/fans, straits, small mountain systems, narrow shelves,
  compact marginal seas, or island-arc spacing must include the parent process
  object/provenance, not just a sharper-looking raster.
- A feature should be called resolution-limited only when its process object
  exists and the higher-resolution run shows the expected morphology; high
  resolution must not be used to excuse missing parent tectonic objects.
- Initial ETOPO-anchored distribution screening is now implemented as
  `aevum.diagnostics.earth_reference`.  It is not yet a direct ETOPO raster fit:
  it uses broad real-Earth envelopes for hypsometry, continental
  shelf/slope/abyss fractions, coastline/ribbon diagnostics, and highland
  coverage.  Direct ETOPO raster sampling should replace these initial
  envelopes before release-quality claims.
- Object presence alone is not enough to claim realistic geomorphology.  A
  generated world should also move toward the real-Earth distribution envelope;
  current 2500-cell Earth-like output is flagged for low land fraction, high
  land ribbon fraction, excessive mean land elevation, too much land above
  2500 m, and slightly high orogen/plateau coverage.
- A follow-up 8000-cell Earth-like release-gate run confirms that the worst
  residual is not merely low resolution.  At 8000 cells the exposed-land ribbon
  fraction rises to `0.677`, narrow necks rise to `22.52` per 1000 land cells,
  largest-coastline complexity rises to `16.48`, and land/continental width
  p50-p90 remains only `2-4` / `3-5` graph steps.  Ocean shelf, abyss, ridge,
  trench, nearshore-superdeep, and far-ocean-shallow fractions remain broadly
  inside the initial screening envelope.  The next calibration target is
  therefore continent/margin geometry and continental hypsometry, not ocean
  bathymetry or renderer resolution.
- The P20 passive-progradation and object-ownership pass restores
  generated-world `interior_basin` visibility in the 2500-cell Earth-like
  regression without weakening the E1 fixture.  It also exposes a separate
  calibration issue: generated continental sediment thickness can saturate near
  the cap, so future basin/platform work should add sediment-routing and
  depocenter contrast instead of classifying basins from a uniform sediment
  blanket.
- P21 adds object-backed supercontinent-breakup seaways, object-seeded
  breakup path extension, continental sediment depocenter contrast, and mature
  suture/accretionary collage relaxation.  The latest 8000-cell Earth-like run
  is back to warning rather than hard failure: land fraction `0.235`, exposed
  land components `13`, ribbon fraction `0.549`, ocean basins `6`, and mean
  land elevation about `1400 m`.  The sediment microbenchmark now separates
  basin/platform/highland/craton sediment means (`3431/1311/869/886 m`).
  Ribbon fraction and coastline complexity remain warning-level, so the next
  Earth-coverage work is upstream continent/margin/rift object geometry and a
  24000-cell medium morphology audit, not renderer changes.  A 72000-cell run
  should be reserved for occasional small-feature deployment tests after
  profiling and medium-resolution validation.
- The first P21 24000-cell medium audit has now been run and it fails on
  medium-resolution supercontinent geometry: land fraction `0.238`, exposed
  land components `11`, largest land component `0.944`, ribbon fraction
  `0.482`, ocean basins `26`, mean land elevation `2102 m`, p95 land elevation
  `5586 m`, and high-land fraction above `2500 m` `0.277`.  This is important:
  ribbon artifacts improve at higher resolution, but the underlying
  continent/rift/seaway object model still leaves one dominant connected
  landmass.
- P22 now has a passing focused fixture for that failure.  Breakup candidates
  carry partition topology metrics and interior rift-axis candidates, so the
  benchmark selects the continent-splitting interior rift instead of peripheral
  weak necks: largest exposed-land component falls to `0.537`, split balance is
  `0.815`, and interior-rift opening reaches `1.000`.
- The follow-up P22 24000-cell audit has now been run and still fails:
  largest exposed-land component improves from `0.944` to `0.856`, but remains
  above the hard `0.82` release threshold; ribbon is `0.487`, land elevation
  mean is `2131 m`, and p95 is `5585 m`.  This anchors P23 as a generated-
  world residual supercontinent problem rather than a fixture-only issue.  The
  release summary now records tectonic object telemetry for future audits so
  the next iteration can separate missing breakup/rift objects from terrain
  non-opening.  A 72000-cell deployment review remains blocked until a 24000
  audit passes the macro-topology gate.
- P23a/P23b added focused coverage for medium rifted continental components,
  continental divergent-boundary rift lifecycle, and boundary-seeded medium
  breakup objects.  In the 24000 generated-world telemetry, rift-system objects
  now exist (`12`), but they are only `3-4` cells each; breakup-seaway count
  remains `1`, breakup area remains about `0.0037`, and largest exposed-land
  component remains `0.852`.  The next Earth-coverage fix is not higher
  resolution; it is broader rift-corridor expression and component-eligibility
  telemetry.
- P23c adds that component telemetry and a multi-corridor rift fixture.  The
  fixture passes, but the 24000 generated world still has only `1` accepted
  breakup object from `10` candidates in the dominant continental component;
  terrain opens only `0.0036` of global area and largest exposed-land component
  remains `0.852`.  The next Earth-coverage fix is therefore object-to-terrain
  seaway effectiveness and exposed-land bridge detection.
- P24a fixes one object-to-terrain failure mode: a breakup axis that is already
  ocean can still seed propagation into adjacent weak land.  The fixture
  passes, but the generated 24000 world has
  `terrain_breakup_seaway_source_reuse=0` and largest exposed-land component
  remains `0.852`, so P24b must measure and repair terrain candidate
  rejection/effectiveness for the current generated world.
- P24b adds that terrain-attempt and stage telemetry.  It showed the accepted
  breakup object was effective (`largest_share 0.945 -> 0.476`) but
  `_regionalize_ocean_floor` later raised the opened shallow seaway back into a
  land bridge.  Object-backed breakup seaway corridors and one ocean-neighbour
  apron are now protected during ocean-floor regionalization.  The 24000 P12
  audit is now `warn` with `0` failed entries: largest exposed-land component
  `0.486`, land fraction `0.237`, component count `12`, ribbon `0.481`, and
  ocean basins `21`.  Remaining Earth-coverage work moves to P25:
  ribbon/coastline complexity, land-fraction undershoot, excessive land
  elevation, and overbroad suture/LIP highland expression.
- P25a implements the first part of that calibration.  Mature suture/LIP and
  stable-craton highlands are now subdued unless supported by active orogenic
  or volcanic process evidence, and `terrain.province` no longer labels every
  suture/LIP land cell as highland.  In the 24000 Earth-like audit, mean land
  elevation drops from about `2131 m` to `1248 m`, p95 from about `5589 m` to
  `3116 m`, land above `2500 m` from about `27.9%` to `9.8%`, and
  suture/LIP/highland province share from about `67%` to `9.3%`.  The audit
  remains `warn`, with `0` hard failures; residual work moves to P25b:
  land-fraction undershoot, ribbon/coastline complexity, the small mean-
  elevation excess, and slightly high abyss fraction.
- The latest resolution-profile preflight
  (`out_profile_900_2500_earth_reference_20260623`) projects about `1.5 min`
  for `8000`, `4.0 min` for `24000`, and `10.1 min` for `72000` cells on the
  current NumPy/SciPy CPU path.  No optional Numba, CuPy, JAX, or Torch
  acceleration package is visible in the current environment, so acceleration
  should begin with measured CPU hotspots before adding optional backends.
  Graph-heavy bottlenecks should be optimized with cached edge reductions or
  optional Numba before considering dense-array GPU backends; CuPy/JAX should
  only be introduced behind an explicit backend boundary with numerical parity
  checks.
  `72000` cells should be used as an occasional small-landform visibility and
  scaling audit, not as a substitute for fixing the failing 24000 macro-
  topology gate.

## Earth Coverage Microbenchmark Suites

These suites complement the R0-R8 refactor microbenchmarks.

### E1. Continental Interior Suite

Initial implementation status, 2026-06-23:

- `terrain.continental_landforms` now derives continental landform objects
  from `terrain.continental_detail`, crust domain/stability/thickness,
  orogen age, sediment thickness, rift potential, and overlapping or adjacent
  tectonic objects.
- Implemented object kinds include `shield`, `platform`, `interior_basin`,
  `rift_basin`, `orogen`, `old_subdued_orogen`, `plateau`,
  `foreland_basin`, and `arc_microcontinent`.
- `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E1 --out out_bench_e1_20260623`
  passes `4/4`.
- Current E1 metrics: basin mean elevation `297 m` vs platform `685 m`;
  basin sediment `1850 m` vs platform `160 m`; recent orogen uplift
  `630 m` vs old orogen uplift `138 m`; collision plateau mean elevation
  `3553 m` vs platform `1285 m`; foreland basin sediment `2400 m` and linked
  parent tectonic object count `2`.
- Remaining limitation: this is fixture-level object coverage.  Global
  Earth-like renders still need E1-E5 integrated release-gate metrics before
  claiming broad real-Earth morphology coverage.

- `E1.craton_platform_basin`
  - fixture: old craton, platform margin, mild extension and sedimentation;
  - metrics: craton stability, platform low relief, interior basin subsidence;
  - acceptance: craton persists, platform remains low relief, basin has
    sediment and lower elevation than adjacent platform.
- `E1.old_orogen_decay`
  - fixture: inactive collision belt after uplift phase;
  - metrics: relief decay through time, preserved suture identity;
  - acceptance: old orogen remains traceable but loses high peak elevations.
- `E1.collision_plateau`
  - fixture: broad thickened collision crust beside lower continental platform;
  - metrics: plateau detail coverage, mean elevation, area, parent object links;
  - acceptance: plateau forms from thickened collision/LIP highland support, not
    by relabeling a one-cell mountain line.
- `E1.foreland_basin`
  - fixture: growing orogen beside stable platform;
  - metrics: subsidence and sediment accumulation on foreland side;
  - acceptance: basin forms adjacent to orogen, not randomly in continent
    interior.

### E2. Margin And Shelf Suite

Initial implementation status, 2026-06-23:

- `terrain.margin_landforms` now derives margin landform objects from
  `ocean.margin_type`, `ocean.depth_province`, `ocean.shelf_width`,
  sediment thickness, river mouth/discharge data, and tectonic boundary or
  lifecycle objects.
- Implemented object kinds include `passive_margin_wedge`, `trench`,
  `forearc_accretionary_prism`, `volcanic_arc`, and `delta_fan`.
- `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E2 --out out_bench_e2_20260623`
  passes `3/3`.
- Current E2 metrics: passive margin shelf p75 `660 m`, slope p50 `1847 m`,
  rise p50 `3186 m`, shelf sediment `1800 m` vs rise sediment `441 m`;
  active trench mean depth `4484 m`, forearc mean elevation `-521 m`, arc
  mean elevation `483 m`; delta/fan sediment `2600 m` and centroid within
  `7.1 deg` of the river mouth.
- Remaining limitation: this is fixture-level object coverage.  Full
  generated-world renders still need E1-E5 integrated release-gate metrics
  before claiming broad real-Earth ocean-basin and arc/plume morphology.

- `E2.passive_margin_shelf_wedge`
  - fixture: rifted margin with sediment supply;
  - metrics: shelf width, slope gradient, rise depth, sediment thickness;
  - acceptance: smooth offshore deepening and sediment wedge growth.
- `E2.active_margin_trench_arc`
  - fixture: ocean-continent subduction;
  - metrics: trench offshore, forearc low, volcanic arc landward, shelf narrow;
  - acceptance: trench and arc remain on the correct sides of polarity.
- `E2.delta_fan`
  - fixture: large river entering passive margin;
  - metrics: delta/fan object area, offshore sediment plume, shelf progradation;
  - acceptance: fan appears at river mouth and not along arbitrary coast.

### E3. Ocean Basin Suite

Initial implementation status, 2026-06-23:

- `terrain.ocean_fabric` now derives ocean-basin fabric objects from
  `ocean.depth_province`, `terrain.elevation_m`, `sediment.thickness_m`,
  `crust.age_myr`, ridge/transform/trench boundary objects, and spreading
  center lifecycle objects.
- Implemented object kinds include `spreading_center`, `transform_fault`,
  `fracture_zone`, `abyssal_plain`, and `age_isochron`.
- `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E3 --out out_bench_e3_20260623`
  passes `3/3`.
- Current E3 metrics: spreading center mean age `4.5 Myr`; transform mean age
  `51.8 Myr`; fracture-zone mean age `107 Myr`; combined fracture-zone span
  `112.7 deg` with `7.5 deg` lat span; abyssal plain mean age `123.6 Myr`,
  sediment `1588 m`, relief p90-p10 `171 m`; age isochron bands are monotonic
  in both age and ridge distance.
- Remaining limitation: this is fixture-level object coverage.  Full
  generated-world renders still need E1-E5 integrated release-gate metrics
  before claiming broad real-Earth ocean-basin, arc/plume, and ice
  geomorphology coverage.

- `E3.ridge_transform_fracture_zone`
  - fixture: offset ridge segments and shear boundary;
  - metrics: transform length, fracture-zone continuity, ridge segment
    continuity, age offsets;
  - acceptance: ocean floor shows continuous ridge-transform geometry, not
    dotted province cells.
- `E3.abyssal_plain_sedimentation`
  - fixture: old ocean basin away from ridge and trench;
  - metrics: depth distribution, sediment smoothing, reduced roughness;
  - acceptance: abyssal plain is broad and smooth but not artificially shallow.
- `E3.ocean_age_isochrons`
  - fixture: stable spreading center through time;
  - metrics: age increases away from ridge and with time, age bands are
    approximately symmetric when spreading is symmetric;
  - acceptance: age derives from spreading history, not only current ridge
    distance.

### E4. Arc, Terrane, Plume Suite

Initial implementation status, 2026-06-23:

- `terrain.arc_plume_landforms` now derives arc/plume landform objects from
  `terrain.elevation_m`, crust origin/domain/thickness, volcanism age,
  active-margin/trench process masks, terrane objects, plume objects, LIP
  objects, and volcano objects.
- Implemented object kinds include `island_arc`, `accreted_terrane`,
  `back_arc_basin`, `hotspot_track`, and `large_igneous_province`.
- `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E4 --out out_bench_e4_20260623`
  passes `4/4`.
- Current E4 metrics: island arc mean elevation `760 m`, accreted terrane
  mean thickness `33000 m` with parent continent id `0`; back-arc basin mean
  age `18 Myr`, mean depth `1441 m`, and landward centroid relative to the
  arc; hotspot distance-age correlation `0.999` with oldest cells much farther
  from the plume than youngest cells; LIP mean crust thickness `43600 m`, mean
  elevation `1500 m`, and parent plume/LIP object count `2`.
- Remaining limitation: this is fixture-level object coverage.  Initial E1-E5
  world-scale release-gate coverage now exists, but generated worlds still
  need stronger parent-object causality and fewer fixture-only feature classes
  before claiming broad real-Earth geomorphology quality.

- `E4.island_arc_accretion`
  - fixture: intra-oceanic arc converges with continent;
  - metrics: arc object, accretion event, parent continent id, suture/terrane
    state;
  - acceptance: terrane attaches only through collision/accretion process.
- `E4.back_arc_basin`
  - fixture: retreating subduction with extensional backarc;
  - metrics: back-arc basin birth, young oceanic crust, arc/backarc geometry;
  - acceptance: basin forms behind arc, not at random ocean cells.
- `E4.hotspot_track`
  - fixture: fixed plume potential under moving plate;
  - metrics: age-progressive volcanic chain, plate-motion alignment, decaying
    activity away from plume head;
  - acceptance: chain is ordered by age and follows plate motion.
- `E4.large_igneous_province`
  - fixture: plume head under continent or near rift;
  - metrics: LIP area, timing, uplift/magmatism, relation to plume potential;
  - acceptance: LIP occurs at plume-potential maximum.

### E5. Ice And Surface Process Suite

Initial implementation status, 2026-06-23:

- `terrain.cryosphere_landforms` now derives cryosphere/surface-process
  objects from `terrain.elevation_m`, `sediment.thickness_m`, current
  `cryosphere.ice_sheet`, optional previous-ice fields, and glacial erosion
  fields when present.
- Implemented object kinds include `ice_sheet_loading`, `glacial_erosion`,
  and `postglacial_rebound`.
- `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite E5 --out out_bench_e5_20260623`
  passes `3/3`.
- Current E5 metrics: ice-sheet loading mean ice thickness `1491 m`, maximum
  ice thickness `2523 m`, estimated bed depression `417 m`; glacial erosion
  mean erosion `960 m`, current ice `900 m`, previous ice `950 m`, relief
  p90-p10 `508 m`; postglacial rebound mean previous ice `1489 m`, unloaded
  ice `1399 m`, and estimated rebound potential `420 m`.
- Remaining limitation: this is fixture-level object coverage.  It does not
  redesign climate, ice dynamics, sea-level response, or seasonal circulation;
  those remain paused until the geomorphology stack has integrated
  release-gate coverage.

- `E5.ice_sheet_loading`
  - fixture: polar continent under cold climate;
  - metrics: ice thickness, bed depression, exposed vs ice-covered terrain;
  - acceptance: ice sheet forms where climate and topography allow.
- `E5.glacial_erosion`
  - fixture: mountain belt under repeated glaciation;
  - metrics: relief modification, valley incision proxy, sediment export;
  - acceptance: glacial modification follows ice flow and relief.
- `E5.postglacial_rebound`
  - fixture: ice sheet decay;
  - metrics: bed rebound, sea-level/coast response;
  - acceptance: rebound follows ice unloading rather than random uplift.

## Earth-Like Release Coverage Gate

Initial implementation status, 2026-06-23:

- P12 release summaries now include an `earth_geomorphology_coverage` section
  per entry and a compact `earth_geomorphology_benchmarks` summary for E1-E5.
- The current integrated Earth-like run at
  `out_p12_earthlike_2500_p18_final_sync_detail_b_20260623/p12_tectonics_release_summary.json`
  reports `24` world-visible features, `3` fixture-only partial
  features, `0` weak features, and `0` missing features.
- Generated parent-object linkage was tightened after the first R8 run:
  `parentless_major_landform_fraction` is now `0.0` in the 2500-cell
  Earth-like coverage run, below the current hard threshold of `0.10`.
- R8 ocean-fabric fallback now promotes generated-world `transform_fault` and
  `fracture_zone` expression in Earth-like runs without explicit transform
  boundary objects.
- R8 continental-object fallback now promotes generated-world
  `interior_basin`, `foreland_basin`, and `old_subdued_orogen` expression from
  stable platform depocenters, recent orogenic loading, and eroded old-orogen
  envelopes.
- P16/R8 plateau classification now promotes generated-world `plateau`
  expression from broad thickened collision/LIP highland cores, with a
  resolution-aware minimum plateau-core width.
- P16/P18 final-step scheduler refresh keeps final terrain/object coverage in
  sync with `t_end`, and process-belt continental detail classification removes
  the previous generated-world orogen/plateau overpaint warning while keeping
  `shield`, `interior_basin`, `orogen`, `old_subdued_orogen`, `plateau`, and
  `transform_fault` world-visible in the 2500-cell Earth-like gate.
- Current warning-level gaps are fixture-only generated-world expression for
  E5 cryosphere landforms, plus the older supercontinent/ribbon morphology
  warnings.  Modern-scale ancient craton calibration remains a separate
  4500 Myr Earth-analogue target, not a current 2500 Myr smoke-gate warning.

The `earth_geomorphology_coverage` section includes at least:

- `covered_feature_count`;
- `missing_feature_count`;
- `weak_feature_count`;
- `parentless_major_landform_fraction`;
- `ridge_transform_continuity_score`;
- `passive_margin_profile_score`;
- `active_margin_profile_score`;
- `ocean_age_isochron_score`;
- `orogen_lifecycle_score`;
- `sedimentary_margin_score`;
- `hotspot_track_score`;
- `ice_geomorphology_score`.

Coverage levels:

- `none`: no object or process exists;
- `weak`: rendered signal exists, but parent object or lifecycle is incomplete;
- `partial`: object and benchmark exist, but release integration is incomplete;
- `world`: object exists in the generated reference world with nontrivial
  area;
- `covered`: future stricter state where object, lifecycle, benchmark,
  diagnostics, and rendered assets all pass.

Initial target:

- no required feature remains `none`;
- all plate-boundary, ocean-basin, and cryosphere/surface-process features
  reach at least `partial`;
- parentless major generated-landform area remains below the current hard
  threshold and continues trending toward zero;
- no major topography in `elevation.png` should ultimately lack a parent
  tectonic or surface-process object.

Long-term target:

- all non-ice required features reach `covered`;
- ice features reach at least `partial` once climate/cryosphere work is active;
- Earth-like maps can be compared to real Earth by feature class, not by visual
  resemblance alone.

## 2026-06-24 P25b Earth-Reference Planform Update

P25b moves the 24000-cell Earth-like audit from a land-budget failure mode
toward a narrower coastline-geometry problem.  Modern coastline smoothing now
pays drowned and unsupported-island cleanup area back into same-component
shallow continental shelf cells, including guarded attached terrane and suture
shelves, while keeping object-backed breakup seaways protected.

Latest 24000-cell audit:

- output:
  `out_p25b_earthlike_24000_planform_payback_final2_20260624/p12_tectonics_release_summary.json`;
- release decision: `warn`, `0` failed, `1` warned;
- land fraction: `0.253`, now inside the initial `0.25-0.33`
  Earth-reference screening range;
- exposed land components: `4`, down from `13` after P25a;
- largest exposed-land component: `0.486`, still modern-fragmented rather than
  single-supercontinent;
- ribbon fraction: `0.410`, improved from about `0.490` after P25a but still
  above the current `0.35` screening target;
- largest-landmass coastline complexity: `21.46`, improved but still far above
  the current `8.0` screening target;
- mean land elevation: `1134 m`, now inside the current broad land-elevation
  screening range, with land p95 `2958 m`;
- remaining Earth-reference misses: ribbon fraction, largest-landmass
  coastline complexity, and abyss fraction.

The important interpretation is that P25b should not be extended into another
global cell-flip cleanup.  The remaining mismatch is evidence that upstream
continental crust, margin, and coastline-width geometry are still too narrow
and jagged.  P26 should therefore work on object-level continent/margin/rift
geometry and ocean-province calibration, then validate with the same 24000
audit and eventual higher-resolution checks.

## 2026-06-24 P26 Spike Rejected

A first P26 attempt widened passive-margin progradation and continental
shape-maintenance candidates.  It passed isolated fixtures, but failed
integration:

- 8000-cell audit initially failed, then still warned with land ribbon
  `0.535`;
- 24000-cell audit failed hard with land `0.216`, components `16`, ribbon
  `0.565`, coastline complexity `34.53`, and `40` ocean basins;
- the change was reverted;
- the restore audit
  `out_p25b_restore_earthlike_24000_after_p26_revert_20260624` returns to the
  accepted P25b state: land `0.253`, components `4`, ribbon `0.410`, coastline
  complexity `21.46`, basins `14`, `warn` with `0` failures.

The lesson for coverage is that P26 needs world-level gates, not only isolated
mechanism fixtures.  Any accepted continent/margin geometry change must keep
the P25b 24000 metrics from regressing while improving ribbon/coastline
complexity and preserving ocean-basin/seam continuity.

This is now enforced by `aevum.diagnostics.p26_regression_gate`.  For any P26
candidate, compare the candidate 24000-cell P12 summary against
`out_p25b_restore_earthlike_24000_after_p26_revert_20260624/p12_tectonics_release_summary.json`.
The gate must pass before the candidate can be counted as progress toward
Earth-reference coverage.

P12 summaries also now include `p26_ribbon_drivers`, a component-level
attribution block for the remaining ribbon/coastline problem.  It lists the top
exposed-land and continental-crust components by ribbon contribution, then
breaks each component down by crust domain, crust origin, continental-detail
province, stability, age, width, and coastline complexity.  This is the
diagnostic entry point for the next P26 mechanism change: first identify
whether the residual mismatch is dominated by young accretionary margins,
suture collages, over-complex stable coastlines, exposed oceanic arc islands,
or narrow articulation necks; then change that upstream process and re-run the
regression gate.

The first attribution-guided 8000-cell audit is
`out_p26_driver_audit_earthlike_8000_20260624`.  It showed that the residual
ribbon/coastline problem is dominated by connected continental crust labelled
as young suture and accreted terrane, not by exposed oceanic crust or by the
hex compiler.  Two narrow mechanism probes were tested and rejected: widening
quiet suture/accretionary maturity worsened ribbon to `0.465`, and adding a
simple active-boundary localization worsened the audit into a hard failure
with ribbon `0.575`.  Both behavior changes were removed; the accepted state is
P25b behavior plus the new read-only attribution and regression-gate tooling.
The attribution block has since been extended with time-since-rework,
recent-orogeny, recent-volcanism, active-rework, and quiet inherited
arc/suture shares so the next P26 change can distinguish currently active
tectonic belts from inherited old collage.
A fresh 8000-cell audit,
`out_p26_time_attribution_earthlike_8000_20260624`, shows that the residual
ribbon problem is currently dominated by recent active rework: the
continental-crust ribbon body has active-rework share about `0.91` and quiet
inherited arc/suture share about `0.00`.  The next Earth-coverage fix should
therefore constrain the footprint of current collision/subduction reworking
rather than simply maturing old collage faster.
`aevum.diagnostics.p26_rework_footprint` and the diagnostic-only `P26`
microbenchmark suite now provide that footprint test: broad recent
collision/suture swaths must be flagged, while localized active belts must be
accepted.  P12 summaries include the footprint metrics for future world-level
audits.
Two active-rework-core production probes were tested after this diagnostic was
added and both were rejected.  One improved visible exposed-land ribbon from
`0.424` to `0.371` but fragmented continental crust and raised continental
ribbon; the other failed hard with ribbon `0.634`.  This means the next
Earth-coverage production change must preserve continental-crust topology and
ancient/stable-crust fractions, not merely narrow the active provenance mask.
The follow-up state layer now exists: `tectonics.deformation_intensity`,
`tectonics.deformation_style`, and `tectonics.deforming_networks` encode active
deforming cores and shoulders separately from crust provenance.  They do not
yet drive terrain or crust-origin changes, but they provide the safer input for
the next Earth-coverage repair.  P12 summaries also now include
`p26_deforming_networks`, which measures active deformation footprint, core /
shoulder split, and overlap with exposed-land and continental ribbon masks.
The deformation state now thins collision, subduction, rift, and transform
axes to plate-contact networks before terrain consumes them; broad
collision/subduction process cells outside the axis remain lower-intensity
shoulders.  The first accepted terrain consumption is semantic only:
`terrain.province`, `terrain.continental_detail`, and terrain landform object
parents read deformation state, while elevation/topographic relief remains on
the P25b path until a separately gated P26 increment can pass 8000/24000
regression.  A first constrained relief increment has now passed the
same-configuration 24000-cell P26 gate: it only raises already exposed, broad,
non-cratonic continental deformation cores/shoulders and leaves land/sea
topology unchanged.  The accepted P26 version additionally gates relief by
exposed-land width and pre-existing elevation (`>=900 m`), with microbenchmark
negative controls proving that active but narrow land cells and active but
low-elevation continental cells receive no relief.  This is a conservative
bridge toward object-level terrain generation, not the final Earth
geomorphology model.
P27 begins that object-level bridge by adding `terrain.orogenic_load` and
`terrain.foreland_accommodation`.  Collision/subduction deformation,
suture/contact processes, and expressed topographic orogens now feed one
shared response state; continental detail and `foreland_basin` landform objects
consume that state and preserve parent tectonic ids.  The attempted
sediment/elevation feedback from this state was rejected because it reduced the
24000-cell Earthlike land budget below the P26 gate, so it remains deferred
until sediment mass balance has its own benchmark.
P28 now supplies that benchmark gate without changing production terrain:
foreland and passive-margin sediment-coupling fixtures must conserve total
continental sediment volume, move sediment into the correct accommodation
zones, draw it from highland/stable source areas, and pass a conservative
projected land-mask check.  This lets the next production step connect sediment
to terrain response with a clear failure signal instead of repeating the P27
land-budget regression.
The next coverage gap after the P28 sediment-budget gate is inland
geomorphology complexity.  Current interiors
can look too flat and monotone even though real continental hypsometry has
broad lowland/platform modes.  P29 should add benchmarked shields, platforms,
interior basins, residual old orogens, rift valleys, plateau margins, volcanic
surfaces, and escarpments before any new visual texture is added.
The first P29 gate now exists as a benchmarked production increment.  It proves
the fixture-level distinction between a mechanism-rich interior, a flat
single-elevation interior, and high-relief checkerboard/speckle texture without
tectonic parentage.  The mechanism-rich fixture also requires parented shield,
old-subdued-orogen, rift-basin, and plateau objects.  Production terrain now
adds a constrained inland-relief pass over already exposed broad continental
interiors, raising old subdued orogens and plateau margins while depressing
interior/rift basins without changing land/ocean topology.  The generated-world
summary also exists in P12 as `p29_inland_geomorphology`.
The first 3000-cell Earthlike smoke reports broad lowland and highland-tail
signals without flagging monotone flat inland, but it is a low-resolution smoke;
the production relief increment is intentionally conservative.  The 8000-cell
follow-up
`out_p29_inland_relief_platform_quantile_earthlike_8000_t2500_20260624`
passes with `0` hard failures and does not regress the accepted planform, but
its inland p25/p50/p75 elevations remain clustered near `1087/1090/1123 m`
and inland IQR remains only `36 m`.  This confirms that broad lowland/platform
hypsometry is Earth-like, but Aevum still over-expresses it as a single smooth
platform.  A wider old-orogen postprocess attempt did not improve the IQR and
regressed ribbon/basin metrics, so the next coverage work should move upstream
into continental-surface regionalization and persistent inland province state,
not broader final-map texture.  The remaining generated Earthlike problems are:
land is still too ribbon-like, the largest coastline remains over-complex,
plateau expression is weak in generated worlds, continental interiors lack
enough shield/platform/basin/old-orogen/rift/escarpment relief partitioning,
and sediment mass-balance terrain feedback remains deferred behind the P28
gate.
P30 now restores the upstream inland-state fixture after separating moderate
platform swells from plateau margins, and P31 adds a dedicated planform
reference gate.  `aevum.diagnostics.planform_reference` can read a generated
world or a P12 release entry and compare land fraction, major landmass
partition, ribbon fraction, and coastline complexity against the initial
real-Earth envelope.  The P31 fixture suite separates a valid broad
multi-continent planform from narrow ribbon land and from broad but jagged
coastline excess.  On the current 24000-cell P29 baseline, only two planform
metrics remain outside envelope: exposed-land ribbon `0.424` versus target max
`0.35`, and largest-landmass coastline complexity `23.276` versus target max
`8.0`.  This makes the next coverage priority explicit: fix coastline
complexity first, then reduce ribbon fraction without collapsing the world into
a single supercontinent or deleting tectonically parented island arcs.
P32 has now added the first coastline-complexity production fixture and
object-constrained local simplification guard.  Its microbenchmarks pass:
production smoothing reduces fixture coastline complexity from `12.63` to
`11.71` and ribbon fraction from `0.363` to `0.329`, while protected breakup
seaways remain ocean.  The current 8000-cell generated world
`out_p32_coastline_complexity_earthlike_8000_t2500_20260624` still warns, but
with no hard failures: land fraction `0.260`, component count `3`, largest land
component `0.648`, ribbon `0.486`, and largest coastline complexity `15.443`.
Nearshore superdeep water (`0.0`) and far-ocean shallow water (`0.089`) are
inside the current screening envelope, so the active blocker has narrowed to
land planform: major-landmass partition, ribbon reduction, coastline
complexity, and generated-world orogen/plateau expression.

## Implementation Dependency

This coverage contract depends on the tectonics refactor plan:

- R1 creates potential fields for plumes, rifts, and thermal evolution;
- R2 gives plate motion physical torque causes;
- R3 creates persistent boundary objects;
- R4 creates ocean-basin and Wilson-cycle lifecycles;
- R5 creates continent and margin lifecycles;
- R6 derives terrain from objects;
- R7 tunes parameters with benchmarks;
- R8 integrates the coverage gate into release validation.

Until those stages exist, the P15 output should be treated as a stable
development fixture, not as proof that Earth geomorphology is fully covered.
