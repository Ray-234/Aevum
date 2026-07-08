# Continental Physiographic Architecture Plan

Status: active research and implementation plan
Owner: tectonics / terrain / diagnostics integration
Created: 2026-06-26
Current entry point: P103 planform mechanism repair

This document archives the next plate-engineering direction after the P69 audit:
the generated worlds can pass many structural and bathymetry guards, but major
continents still tend to express as broad, high, smooth platforms.  The missing
layer is not another final terrain cleanup.  The engine needs a process-derived
continental physiographic architecture: large continents should be assembled
from multiple tectonic and geomorphic provinces, and those province objects
should control elevation, sediment, drainage, and final landform compilation.

The expanded research, source-collection, theory, optimization, and
microbenchmark plan for this work is archived in
`docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`.
The 2026-06-27 evidence-collection expansion in that document is the canonical
archive for staged source packets, theory/test details, implementation methods,
and optimization targets.  The P101+ and P103+ real-Earth comparison archives in
the same document are the current execution plan for clearing the remaining
planform and crust/sediment residual blockers without weakening the reference
envelopes.

## Problem Statement

The current terrain pipeline builds most exposed continental relief from crustal
isostasy, broad continental support, and later smoothing.  `terrain.province`,
`terrain.continental_detail`, and continental landform objects then classify or
lightly adjust a surface that has already been made.  This ordering can pass
coarse land-elevation envelopes while still producing unrealistic interiors:

- large inland areas sit at similar elevation;
- platforms are over-expressed as a single smooth highland or tableland;
- internal basins, old orogens, rift systems, foreland basins, shields, and
  passive-margin lowlands do not reliably partition each major continent;
- province boundaries exist mainly as raster field thresholds, not persistent
  region objects with parent processes and adjacency;
- current microbenchmarks prove local recognition and response, but do not yet
  prove that generated continents contain Earth-like internal province graphs.

P69 evidence for the Earth-like reference member:

- land fraction `0.248`;
- land component count `2`;
- largest land component fraction `0.614`;
- land ribbon fraction `0.412`;
- largest coastline complexity `17.44`;
- abyss fraction `0.687`;
- P48 mature cratonic share `0.254`, with basin share `0.046`.

Those values show that the modern Earth-like blocker has shifted from local
bathymetry bugs toward continent architecture and real-Earth calibration.

## Target Model

Each major exposed continent should become a graph of process-derived
physiographic provinces before final terrain and map compilation.

Core province types:

- `shield`: exposed ancient cratonic basement, old stable lithosphere, subdued
  relief, resistant basement texture.
- `covered_platform`: old cratonic basement below sediment cover, broad low to
  moderate relief, common continental interior mode.
- `intracratonic_basin`: long-lived sag or warping on a mature platform,
  low elevation, high sediment, usually not tied to active margins.
- `old_orogen`: old collision or accretion belt, inherited suture, eroded but
  still a natural boundary and moderate relief source.
- `active_orogen`: active or recent convergent or collisional belt, high relief,
  narrow to moderate width, strong parent boundary.
- `orogenic_plateau`: thickened crust behind or inside active/recent collision,
  high but area-limited.
- `foreland_basin`: flexural accommodation adjacent to an orogenic load, low
  elevation and high sediment.
- `rift_basin`: extensional axis, linear lowland or lake-chain candidate.
- `rift_shoulder`: uplifted shoulder adjacent to a rift basin.
- `volcanic_lip_province`: plume/LIP/large volcanic surface, local plateau or
  swell, finite age and decay.
- `passive_margin_lowland`: mature rifted margin, coastal plain, shelf-linked
  sediment wedge and low relief.
- `escarpment`: boundary between platform, rift shoulder, plateau, or coastal
  plain, produced by differential uplift and erosion.

Boundary sources:

- old sutures and mobile belts;
- craton and shield margins;
- inherited rift axes and failed rifts;
- active plate boundaries and deformation zones;
- foreland flexural hinge lines;
- passive-margin hinge zones and coastal plains;
- lithosphere thickness / thermal support contrasts;
- sedimentary basin margins;
- drainage divides and long-term erosional escarpments.

## Research Corpus

The research corpus should be collected in stages.  Prefer sources that can be
turned into reproducible metrics or fixture definitions.

### A. Topography, Bathymetry, and Planform

- NOAA ETOPO 2022 Global Relief Model:
  https://www.ncei.noaa.gov/products/etopo-global-relief-model
  - Use for land hypsometry, ocean hypsometry, shelf/slope/abyss references,
    highland tail, lowland share, and global rendered comparisons.
- GEBCO Gridded Bathymetry:
  https://www.gebco.net/data-products/gridded-bathymetry-data
  - Use for bathymetry cross-checks, ocean basin/shelf/slope classification,
    and ridge/trench/deep basin proportions.
- Natural Earth physical vectors:
  https://www.naturalearthdata.com/downloads/10m-physical-vectors/
  - Use for land polygons, coastline complexity, island/component counts,
    rivers/lakes at cartographic scale, and map-compiler QA.

### B. Physiographic and Geological Provinces

- USGS Physiographic Divisions of the conterminous United States:
  https://data.usgs.gov/datacatalog/data/USGS%3Ae04ea9e9-17b6-45ae-b279-7bc35ea79539
  - Use as a high-quality labeled case study: divisions, provinces, and
    sections share topography, rock type, structure, and geologic history.
- NPS physiographic province explanation:
  https://www.nps.gov/subjects/geology/physiographic-provinces.htm
  - Use for acceptance language: provinces are defined by geomorphology,
    geologic structures, climate, underlying geology, and geologic history.
- Global geological province and tectonic plate compilation:
  https://zenodo.org/records/6586972
  - Use for modern global province classes, active plate boundaries,
    deformation zones, orogens, shields, cratons, and province labels.
- USGS world petroleum/geologic provinces:
  https://data.usgs.gov/datacatalog/
  - Use cautiously as supporting geologic province boundaries, especially for
    basins; some boundaries are assessment-oriented rather than strictly
    physiographic.

### C. Plate Reconstructions and Tectonic Forcing

- GPlates:
  https://www.gplates.org/
  - Use for deep-time plate reconstruction workflows and visual validation.
- pyGPlates / EarthByte:
  https://www.earthbyte.org/category/resources/software-workflows/pygplates/
  - Use for future automated interrogation of plate reconstructions and
    reconstructable feature geometries.
- EarthByte plate motion and reconstruction data:
  https://www.earthbyte.org/category/reconstructions/
  - Use for Wilson-cycle case studies and reference time slices.
- PB2002 plate boundary model:
  https://peterbird.name/oldFTP/PB2002/
  - Use for present-day plate boundary types and rough process localization.
- GEM Global Active Faults:
  https://github.com/GEMScienceTools/gem-global-active-faults
  - Use for modern active fault density and deformation belt references.

### D. Crust, Lithology, Sediment, and Drainage

- CRUST1.0:
  https://ds.iris.edu/ds/products/emc-crust10/
  - Use for crustal thickness, sediment thickness, and broad continental /
    oceanic crust checks.
- GLiM Global Lithological Map:
  https://www.geo.uni-hamburg.de/en/geologie/forschung/aquatische-geochemie/glim.html
  - Use for lithology class references and surface expression of shields,
    platforms, volcanic provinces, carbonate platforms, and sedimentary basins.
- NOAA total sediment thickness:
  https://www.ncei.noaa.gov/products/total-sediment-thickness-oceans-seas
  - Use for ocean and marginal sea sediment-thickness calibration.
- HydroSHEDS / HydroBASINS / HydroRIVERS:
  https://www.hydrosheds.org/products
  - Use for drainage basins, river networks, flow accumulation, and watershed
    shape metrics.
- GMBA mountain inventory:
  https://www.earthenv.org/mountains
  - Use for mountain range count, size, hierarchy, and regional mountain
    distribution checks.

## Theory Work Packages

### T1. Continental Assembly and Wilson Cycle

Questions:

- Which province objects should form during rifting, passive-margin maturity,
  subduction initiation, ocean closure, collision, and post-orogenic erosion?
- Which objects persist after the active boundary is gone?
- Which objects may be reactivated during later breakup or collision?

Required outputs:

- lifecycle table for each province type;
- parent-process requirements;
- age-decay rules;
- accepted and forbidden transitions.

### T2. Craton, Shield, Platform, and Basin Logic

Questions:

- How much of a mature continent should be old stable basement versus exposed
  shield versus sediment-covered platform?
- When should a platform become an intracratonic basin?
- How should stable cratonic lithosphere affect elevation without making the
  entire interior high and flat?

Required outputs:

- shield/platform/basin object-generation rules;
- lowland and highland envelope by province;
- sediment coupling rules;
- erosion and peneplain representation.

### T3. Orogens, Plateaus, Forelands, and Old Sutures

Questions:

- How should active orogens differ from old orogens in relief, width, and
  persistence?
- How should foreland basins be placed relative to mountain loads and cratons?
- When can thick crust produce a plateau rather than a narrow mountain belt?

Required outputs:

- active-orogen and old-orogen relief kernels;
- foreland-basin adjacency constraints;
- plateau area caps and decay rules;
- suture boundary preservation.

### T4. Rift and Passive-Margin Provinces

Questions:

- How should a rift axis, shoulder, failed rift, and passive-margin lowland
  evolve through time?
- Which rifts become ocean gateways and which remain internal basins?
- How do rifted margins create shelves, coastal plains, and sediment wedges?

Required outputs:

- rift-axis and shoulder object rules;
- failed-rift/aulacogen preservation;
- passive-margin lowland and shelf coupling;
- no-random-seaway constraints.

### T5. Drainage, Erosion, and Sediment Feedback

Questions:

- How should drainage divides follow province boundaries?
- How strongly can long-term erosion lower old orogens and cratons?
- How should sediment move from highland sources to foreland, intracratonic,
  passive-margin, and ocean basins?

Required outputs:

- drainage-before-final-detail workflow;
- sediment source-to-sink rules;
- erosion rates by province and climate placeholder;
- basin accommodation accounting.

## Implementation Architecture

The next design should add a first-class province architecture layer between
crust/tectonics and final terrain.

Proposed fields:

- `tectonics.continental_province_code`
- `tectonics.continental_province_id`
- `tectonics.province_boundary_strength`
- `tectonics.province_parent_process`
- `terrain.physiographic_province`
- `terrain.province_relief_template_m`
- `terrain.province_base_elevation_m`
- `terrain.province_erosion_susceptibility`

Proposed object sets:

- `tectonics.shields`
- `tectonics.platforms`
- `tectonics.interior_basins`
- `tectonics.old_orogens`
- `tectonics.active_orogens`
- `tectonics.foreland_basins`
- `tectonics.rift_basins`
- `tectonics.rift_shoulders`
- `tectonics.passive_margin_lowlands`
- `tectonics.volcanic_lip_provinces`
- `terrain.escarpments`
- `terrain.drainage_divides`

Proposed generation order:

1. Generate plate, crust, domain, and persistent basement state.
2. Identify major continental components and stable continent IDs.
3. Build a province graph per major continent from parent processes.
4. Assign province codes, IDs, boundaries, and parent-object links.
5. Generate regional base elevation from province templates.
6. Apply orogen/rift/foreland/passive-margin relief kernels.
7. Route drainage and sediment using province boundaries and base relief.
8. Apply erosion and sediment accommodation.
9. Compile `terrain.continental_detail` and landform objects from the province
   graph, not directly from final elevation thresholds.
10. Run map compiler and release diagnostics.

Hard design constraints:

- A major continent must not be represented by one dominant platform object
  unless it is intentionally a small, simple continent.
- Active orogens and plateaus must be area-limited and parented.
- Basins and lowlands must have explicit accommodation causes.
- Old sutures should persist as boundaries even after relief decays.
- Random seeds may perturb textures, but may not choose the existence or
  location of major province objects.

## Microbenchmark Plan

### P70. Research Corpus and Metric Scaffold

Purpose:

- archive source inventory;
- define real-Earth metric extraction targets;
- add diagnostics that can consume reference data later without requiring the
  data to be bundled into the repo.

Acceptance:

- source inventory is documented;
- metric schema exists for land hypsometry, mountain share, basin share,
  province count, province adjacency, drainage shape, and shelf/slope/abyss;
- current P69 generated metrics are mapped to the new schema.

### P71. Province Graph Fixture Suite

Fixtures:

- `craton_platform_basin_fixture`;
- `active_orogen_foreland_fixture`;
- `old_suture_orogen_fixture`;
- `rift_axis_shoulder_fixture`;
- `passive_margin_lowland_fixture`;
- `volcanic_lip_plateau_fixture`;
- `multi_province_continent_fixture`.

Acceptance:

- each fixture creates the expected province objects;
- parent-process links are present;
- province IDs are spatially contiguous where expected;
- no province is created from random texture alone;
- each fixture produces expected ordering of elevation and sediment.

### P72. Generated-World Province Diversity Gate

Metrics:

- major continent count;
- province count per major continent;
- largest province fraction per major continent;
- basin/platform/shield/old-orogen/rift/passive-margin shares;
- highland share by province;
- inland p90-p10 and p95-p05 relief;
- lowland fraction below 500 m and below 1000 m;
- parented highland fraction.

Acceptance:

- every major Earth-like continent has at least three province classes;
- largest internal province fraction is capped for large continents;
- basin or passive-margin lowland share has a minimum;
- active highland/plateau share has a maximum;
- highlands without tectonic parents remain below a small tolerance.

### P73. Real-Earth Case-Study Calibration

Case studies:

- North America: Canadian Shield, Interior Plains, Appalachians, Cordillera,
  Basin and Range, Coastal Plain.
- South America: Andes, Amazon foreland/basin, Brazilian Shield, Patagonia.
- Africa: West/Central/Southern cratons, Congo Basin, Sahara platform, East
  African Rift, Ethiopian Highlands.
- Eurasia: East European Platform, West Siberian Basin, Alps-Himalaya belt,
  Tibetan Plateau, Central Asian old orogens.
- Australia: Western Shield, Central Lowlands, Eastern Highlands, passive
  margins.

Acceptance:

- each case has a reference province sketch and metric envelope;
- generated worlds can reproduce the feature class, not the exact geography;
- failures are categorized as missing process, wrong amplitude, wrong scale,
  wrong adjacency, or compiler/rendering mismatch.

### P74. Terrain Coupling Rewrite

Purpose:

- move regional base elevation and relief response from late cell-level passes
  into province-object templates.

Acceptance:

- generated high inland plateaus are reduced unless parented by plateau/orogen
  processes;
- old platforms and basins provide broad lowland modes;
- old orogens remain visible as subdued belts;
- rift axes and forelands form lows in the correct adjacency;
- existing P69 bathymetry/compiler guards do not regress.

### P75. Release and Promotion Audit

Purpose:

- rerun a high-resolution ensemble and decide whether P68/P70+ mature continent
  behavior can be promoted toward default Earth-like generation.

Acceptance:

- P69-style assets complete at high resolution;
- Earth-like reference enters the updated real-Earth envelope for planform and
  province architecture;
- wet/dry physical variants stay plausible without being forced to match modern
  Earth exactly;
- old P29/P48/P49/P68 regression gates still pass or are replaced by stronger
  province-architecture equivalents.

## Optimization Targets

Targets should be treated as initial envelopes and revised after reference data
extraction.

Generated Earth-like reference:

- exposed land fraction: close to modern Earth envelope, approximately
  `0.25-0.33`;
- major land components: multiple, not one dominant supercontinent;
- largest land component fraction: below the initial P69 blocker;
- land ribbon fraction: below current P69 reference and wet member values;
- coastline complexity: lower than the P69 reference largest-component value;
- land mean elevation: keep inside broad Earth-like envelope, but avoid using a
  uniform high platform to achieve it;
- highland area above 2500 m: present but area-limited;
- province count per major continent: at least three meaningful classes;
- largest internal province fraction: capped for large continents;
- basin/passive-margin lowland share: positive and visible;
- old-orogen/rift/shield/platform/foreland expression: parented by objects;
- abyss/shelf/slope fractions: preserve P69 bathymetry sanity.

## Testing and Data Policy

- Do not commit large reference rasters into the repo by default.
- Add download or extraction scripts only when licensing permits and the source
  is stable.
- Store derived small metric JSON fixtures when possible.
- Each reference-derived metric file should record source URL, source version,
  extraction date, projection/resolution, and preprocessing steps.
- Benchmarks must run without network by default; network/data acquisition
  should be a separate preparation step.
- Visual contact sheets are diagnostic evidence, not acceptance by themselves.

## Progress Log

2026-06-26 - Plan archived

- Created this plan after the P69 high-resolution audit and the inland-high-flat
  root-cause review.
- Reframed P70 from local planform cleanup to research-corpus and
  province-architecture metric scaffolding.
- Current implementation status:
  - P29 provides conservative inland geomorphology response.
  - P48/P49 provide shield/platform/interior-basin object and subsidence
    coupling fixtures.
  - P68/P69 provide high-resolution physical-ensemble evidence.
  - Missing: first-class multi-province graph per generated continent.

Next:

- P71 should implement deterministic province-graph fixture objects before
  rewriting production terrain.

2026-06-26 - P70 reference scaffold implemented

- Added an offline P70 source inventory and metric schema in
  `aevum/diagnostics/physiographic_reference.py`.
- Mapped the P69 Earth-like high-resolution baseline into the new schema,
  including the six out-of-envelope blockers and the explicit architectural
  gaps: no first-class multi-province graph, no required multi-class province
  diversity per major continent, basin lowland expression not guaranteed, and
  province boundaries not yet primary terrain drivers.
- Added `P70.reference_corpus_metric_scaffold` to the tectonics benchmark CLI
  and pytest suite.  The benchmark deliberately records that real reference
  data download/extraction and production province graph generation are still
  pending.

Next:

- P72: move from fixtures to generated-world gates, requiring every major
  continent to contain multiple process-parented physiographic provinces.

2026-06-26 - P71 province graph fixture suite implemented

- Added deterministic raster-graph fixtures in
  `aevum/diagnostics/physiographic_provinces.py`.
- Covered the planned fixture families:
  `craton_platform_basin_fixture`, `active_orogen_foreland_fixture`,
  `old_suture_orogen_fixture`, `rift_axis_shoulder_fixture`,
  `passive_margin_lowland_fixture`, `volcanic_lip_plateau_fixture`, and
  `multi_province_continent_fixture`.
- Each fixture validates expected province classes, class adjacencies,
  process-parent links, spatial contiguity of province IDs, non-random source
  attribution, and expected elevation/sediment ordering.
- Added `P71.province_graph_fixture_suite` to the tectonics benchmark CLI and
  pytest suite.  The benchmark intentionally records that production generated
  world integration is still pending.

Next:

- P73: add real-Earth case-study calibration sketches now that generated-world
  province diversity has a stable metric surface.
- P74: rewrite production terrain coupling so passive-margin lowlands and
  province boundaries become first-class terrain drivers instead of proxies.

2026-06-26 - P72 generated-world province diversity gate implemented

- Added `aevum/diagnostics/province_diversity.py`, which analyzes a generated
  world by exposed continental component.
- The P72 diagnostic measures, per major component: province/detail-class
  diversity, largest internal province fraction, basin-or-lowland proxy share,
  lowland fractions below 500 m and 1000 m, active highland/plateau share,
  relief, highland parentage proxy, and overlapping terrain landform object
  context.
- Added `P72.generated_world_province_diversity_gate` to the tectonics
  benchmark CLI and pytest suite.  The benchmark runs a 900-cell Earth-like
  world to keep the gate executable while still using the real generator.
- Current generated-world status: the gate passes with multiple major exposed
  continental components, at least three province classes per major component,
  capped largest internal province fraction, parented highlands, and parented
  terrain landform objects.
- Explicit limitation retained for P74: `terrain.continental_detail` still does
  not have a dedicated passive-margin-lowland class, so P72 uses basin/rift
  basin/lowland fractions as a proxy.

Next:

- P75: rerun a release/promotion audit after P74 proves terrain and compiler
  behavior.

2026-06-26 - P73 real-Earth case-study calibration implemented

- Added `aevum/diagnostics/earth_case_studies.py`, an offline case-study
  calibration registry for North America, South America, Africa, Eurasia, and
  Australia.
- Each case now records a reference province sketch, province classes, parent
  processes, expected adjacency edges, broad metric envelopes, source IDs, and
  failure categories.
- The calibration explicitly requires generated worlds to reproduce feature
  classes and relationships, not exact modern geography.
- Added `P73.real_earth_case_study_calibration` to the tectonics benchmark CLI
  and pytest suite.
- Explicit limitations retained for later work:
  - reference raster extraction remains pending;
  - P75 release/promotion audit remains pending.

Next:

- P75: rerun a release/promotion audit after P74 proves terrain and compiler
  behavior.

2026-06-26 - P74 terrain coupling rewrite implemented

- Added a production P74 province-template terrain response in
  `aevum/modules/terrain.py`.
- Added `terrain.passive_margin_lowland`, a first-class process-template field
  for passive-margin coastal plains and lowlands.
- Added `passive_margin_lowland` to `terrain.continental_landforms`, with
  parent process `passive_margin_subsidence_and_coastal_plain_sedimentation`.
- The P74 terrain response now gives process-parented low-amplitude terrain
  adjustments to passive-margin lowlands, foreland lows, rift lows, subdued old
  orogen templates, and unsupported highlands while preserving land/ocean sign.
- Added a conservative Earth-like excess-oceanic-land trim that only drowns
  low exposed oceanic remnants when their share of total exposed land exceeds
  the regression envelope.
- Updated generated-world province diversity so P72 now sees the P74
  passive-margin lowland object layer as active.  The detail enum is still not
  expanded; passive-margin lowland is first-class as a field/object template.
- Added `P74.terrain_coupling_rewrite` to the tectonics benchmark CLI and
  pytest suite.
- Current 900-cell Earth-like P74 gate result:
  - template response area fraction: `0.0389`;
  - passive-margin lowland field area fraction: `0.0078`;
  - passive-margin lowland object count: `2`;
  - terrain landform kind count: `8`;
  - nearshore superdeep and compiler sign mismatch fractions: `0.0`.

Next:

- P75: rerun the release/promotion audit, checking P69-style assets, P70-P74
  gates, bathymetry/compiler stability, and whether the Earth-like reference is
  ready for promotion.

2026-06-26 - Detailed real-Earth geomorphology research plan archived

- Added `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`.
- The new archive expands the P70-P75 direction into staged source collection,
  theory modules, metric extraction, generated-world gates, production
  province-architecture work, drainage/sediment coupling, and high-resolution
  promotion audits.
- It also proposes follow-up suites P76-P89 so future work can enrich tests and
  implementation one process at a time rather than relying on final rendered
  maps.

Next:

- Start Stage A of the archived plan with source-ledger and real-Earth
  hypsometry/province metric extraction.

2026-06-27 - P75 release and promotion audit implemented

- Added `P75.release_and_promotion_audit` and `run_p75_bench` to the tectonics
  benchmark CLI.
- P75 reads the archived P69 high-resolution physical-ensemble audit assets and
  verifies that the selected 8000-cell P69 evidence has complete key PNG asset
  sets.
- P75 also reruns P70-P74 gates and the legacy P29/P48/P49/P68 regression
  gates so the promotion decision is tied to current code, not only archived
  intent.
- Current P75 output:
  `/Users/rayw/Projects/aevum/out_bench_p75_release_promotion_audit_20260627/`.
- Current P75 metrics:
  - selected P69 summary:
    `/Users/rayw/Projects/aevum/out_bench_p69_highres_physical_ensemble_visual_audit_20260626_c/tectonics_bench_summary.json`;
  - P69 cells: `8000`;
  - P69 ready members: `3/3`;
  - P69 key asset files present: `18/18`;
  - P70-P74 gate pass count: `5/5`;
  - P29/P48/P49/P68 legacy gate pass count: `4/4`;
  - promotion ready: `False`;
  - explicit blockers:
    `reference_data_download_or_raster_extraction_pending`,
    `first_class_production_province_graph_pending`, and
    `p69_earthlike_reference_needs_calibration`.
- Verification:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p75_release_and_promotion_audit_pass -q`
  -> `1 passed in 224.72s`.
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P75 --out out_bench_p75_release_promotion_audit_20260627`
  -> `status: pass`.

Next:

- Stage A starts with `P76.reference_source_ledger_schema` and
  `P77.real_earth_hypsometry_extraction`.

2026-06-27 - P76 reference source ledger schema implemented

- Added `aevum/diagnostics/reference_source_ledger.py`.
- Added `P76.reference_source_ledger_schema` and `run_p76_bench` to the
  tectonics benchmark CLI.
- P76 converts the P70 source inventory into an explicit offline ledger schema
  with required fields for source id, category, phase ids, URL, source kind,
  version note, license status, acquisition status, extraction status, raw data
  policy, local storage policy, projection/resolution note, checksum status,
  and derived metric targets.
- Current P76 output:
  `/Users/rayw/Projects/aevum/out_bench_p76_reference_source_ledger_schema_20260627/`.
- Current P76 metrics:
  - source count: `18`;
  - external source count: `17`;
  - internal source count: `1`;
  - category count: `12`;
  - required field count: `14`;
  - missing required field source count: `0`;
  - invalid phase reference source count: `0`;
  - duplicate source id count: `0`;
  - external raw data required by benchmark: `False`.
- Verification:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p76_reference_source_ledger_schema_pass -q`
  -> `1 passed in 0.65s`.
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P76 --out out_bench_p76_reference_source_ledger_schema_20260627`
  -> `status: pass`.

Next:

- Implement `P77.real_earth_hypsometry_extraction` with reproducible small
  derived metric fixtures, still without committing raw ETOPO/GEBCO rasters.

2026-06-27 - P77 real-Earth hypsometry extraction fixture implemented

- Added `data/reference/earth_hypsometry_fixture_20260627.json`, a small
  derived metric fixture for global land/ocean hypsometry and bathymetry
  partition checks.
- Added `aevum/diagnostics/real_earth_hypsometry.py`.
- Added `P77.real_earth_hypsometry_extraction` and `run_p77_bench` to the
  tectonics benchmark CLI.
- P77 validates that the fixture references ETOPO/GEBCO source IDs, does not
  store raw rasters in the repo, has area-weighted bins that sum to one, has
  land/ocean bin totals matching derived metrics, and passes the current broad
  Earth-reference envelopes.
- Current P77 output:
  `/Users/rayw/Projects/aevum/out_bench_p77_real_earth_hypsometry_extraction_20260627/`.
- Current P77 metrics:
  - bin count: `12`;
  - land fraction: `0.292`;
  - land elevation mean: `835 m`;
  - land elevation p95: `2600 m`;
  - high land fraction above 2500 m: `0.055`;
  - shelf/slope-rise/abyss fractions of ocean: `0.056/0.117/0.636`;
  - shelf-to-abyss depth delta: `3900 m`;
  - envelope checks passing: `8/8`;
  - raw raster stored in repo: `False`;
  - direct raster extraction still marked pending: `True`.
- Verification:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p77_real_earth_hypsometry_extraction_pass -q`
  -> `1 passed in 0.66s`.
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P77 --out out_bench_p77_real_earth_hypsometry_extraction_20260627`
  -> `status: pass`.

Next:

- Implement `P78.generated_hypsometry_envelope` so generated worlds are checked
  against the P77 fixture schema rather than only the older screening ranges.

2026-06-27 - P78 generated hypsometry envelope implemented

- Added generated-world hypsometry extraction and fixture comparison helpers to
  `aevum/diagnostics/real_earth_hypsometry.py`.
- Added `P78.generated_hypsometry_envelope` and `run_p78_bench` to the
  tectonics benchmark CLI.
- P78 compares a current 900-cell Earth-like generated world against the P77
  fixture schema and also checks that archived P69 8000-cell evidence remains
  available.
- Current P78 output:
  `/Users/rayw/Projects/aevum/out_bench_p78_generated_hypsometry_envelope_20260627/`.
- Current 900-cell generated-world metrics:
  - land fraction: `0.240`;
  - land elevation mean: `859 m`;
  - land elevation p95: `1928 m`;
  - high land fraction above 2500 m: `0.009`;
  - lowland fractions below 500/1000 m: `0.389/0.532`;
  - shelf/slope-rise/abyss ocean fractions: `0.091/0.190/0.308`;
  - shelf-to-abyss depth delta: `3387 m`;
  - core hypsometry envelope pass: `True`.
- P78 still records expected residual blockers:
  `land_fraction`, `land_component_count`, `land_ribbon_fraction_gt_0_5`,
  `land_coastline_complexity_largest`, and `trench_fraction_of_ocean`.
  Promotion calibration remains required.
- Verification:
  `.venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p78_generated_hypsometry_envelope_pass -q`
  -> `1 passed in 30.85s`.
  `.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P78 --out out_bench_p78_generated_hypsometry_envelope_20260627`
  -> `status: pass`.

Next:

- Implement `P79.province_reference_graph_extraction`.

2026-06-27 - P79 province reference graph extraction implemented

- Added `aevum/diagnostics/province_reference_graph.py`.
- Added `P79.province_reference_graph_extraction` and `run_p79_bench` to the
  tectonics benchmark CLI.
- P79 converts the P73 real-Earth case-study sketches into a compact derived
  province reference graph.  It records province nodes, adjacency edges,
  province classes, parent processes, class-edge coverage, source IDs, and
  extraction policy.
- The fixture is intentionally not raw GIS data.  Raw vector extraction remains
  pending and explicitly marked as such; the current fixture is a small
  executable reference for generated province graphs.
- Current P79 output:
  `/Users/rayw/Projects/aevum/out_bench_p79_province_reference_graph_extraction_20260627/`.
- Current P79 metrics:
  - case count: `5`;
  - province node count: `29`;
  - province adjacency edge count: `25`;
  - province class count: `9`;
  - parent process count: `15`;
  - source ID count: `9`;
  - class-edge count: `16`;
  - missing required feature classes/processes/class edges: `0/0/0`;
  - raw vectors stored in repo: `False`;
  - direct vector extraction still marked pending: `True`.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p79_province_reference_graph_extraction_pass -q`
  -> `1 passed in 0.29s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P79 --out out_bench_p79_province_reference_graph_extraction_20260627`
  -> `status: pass`.

Next:

- Implement `P80.generated_province_graph_reference_comparison`.

2026-06-27 - P80 generated province graph reference comparison implemented

- Added `aevum/diagnostics/generated_province_reference.py`.
- Added `P80.generated_province_graph_reference_comparison` and
  `run_p80_bench` to the tectonics benchmark CLI.
- P80 maps the current generated world's `terrain.continental_detail` and
  `terrain.continental_landforms` into the P79 reference province classes, then
  compares global and major-continent class coverage, parent-process coverage,
  class-edge coverage, dominant-class fraction, and recorded residual blockers.
- Current P80 output:
  `/Users/rayw/Projects/aevum/out_bench_p80_generated_province_graph_reference_comparison_20260627/`.
- Current P80 metrics:
  - cells: `900`;
  - major continent count: `2`;
  - generated/reference class count: `8/9`;
  - mapped/required parent process count: `8/8`;
  - generated/required class-edge count: `23/9`;
  - missing required feature classes/processes: `0/0`;
  - missing required class edges: `1`;
  - unexpected missing reference classes/edges: `0/0`;
  - min major-continent reference class count: `7`;
  - max largest reference-class fraction: `0.473`;
  - landform object/overlay cell count: `31/197`.
- Recorded residuals:
  - missing reference class: `volcanic_lip_plateau`;
  - missing required class edge: `rift_system|volcanic_lip_plateau`.
- Interpretation:
  - The generated world now compares cleanly against P79 for required feature
    classes, required parent processes, and most class adjacency.  Promotion is
    still blocked by the lack of first-class production province IDs and missing
    volcanic LIP/plateau expression.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass -q`
  -> `1 passed in 31.00s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P80 --out out_bench_p80_generated_province_graph_reference_comparison_20260627`
  -> `status: pass`.

Next:

- Implement `P81.boundary_process_geometry_reference`.

2026-06-27 - P81 boundary process geometry reference implemented

- Added `aevum/diagnostics/boundary_process_geometry.py`.
- Added `P81.boundary_process_geometry_reference` and `run_p81_bench` to the
  tectonics benchmark CLI.
- P81 defines a compact boundary-process geometry fixture using PB2002/GPlates/
  EarthByte/GEM source IDs, broad reference length-fraction envelopes, and a
  deterministic spherical synthetic network with ridge, transform, subduction
  trench, collision/suture, diffuse deformation, passive margin, and
  continental rift process types.
- The synthetic fixture checks length shares, transform-offset adjacency to
  ridges, trench-active-margin adjacency, collision-diffuse adjacency, curved
  trench geometry, and antimeridian continuity on the sphere.
- P81 also compares the current generated Earth-like boundary masks and records
  the expected current residual: transform boundaries are still absent.
- Current P81 output:
  `/Users/rayw/Projects/aevum/out_bench_p81_boundary_process_geometry_reference_20260627/`.
- Current P81 metrics:
  - source ID count: `5`;
  - synthetic process type count: `7`;
  - synthetic length-envelope checks passing: `7/7`;
  - synthetic transform offset count: `2`;
  - synthetic antimeridian ridge component count: `1`;
  - synthetic transform-near-ridge fraction: `0.846`;
  - synthetic trench-near-active-margin fraction: `1.000`;
  - current generated process type count: `6`;
  - current generated boundary cells: `103`;
  - current generated boundary cell fraction: `0.114`;
  - current generated missing process types: `transform`.
- Interpretation:
  - The boundary geometry reference and synthetic process-network gate are now
    executable.  The current Earth-like world has ridge, trench/subduction,
    collision/suture, diffuse deformation, passive margin, and continental rift
    masks, but still lacks first-class transform boundary expression.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p81_boundary_process_geometry_reference_pass -q`
  -> `1 passed in 31.65s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P81 --out out_bench_p81_boundary_process_geometry_reference_20260627`
  -> `status: pass`.

Next:

- Implement `P82.wilson_cycle_lifecycle_reference`.

2026-06-27 - P82 Wilson-cycle lifecycle reference implemented

- Added `aevum/diagnostics/wilson_cycle_lifecycle.py`.
- Added `P82.wilson_cycle_lifecycle_reference` and `run_p82_bench` to the
  tectonics benchmark CLI.
- P82 defines a deterministic scripted Wilson-cycle reference covering seven
  stages: continental rift, spreading ocean, mature/passive-margin ocean,
  subduction closure, arc-collision margin, suture relict, and old-orogen
  relict.
- The scripted fixture verifies that the same basin id and lineage key persist
  through all stages, phase codes progress monotonically, basin age increases,
  rift/passive/spreading/closing/suture/Wilson/gateway object sets all appear,
  gateway parent causality is intact, and the old orogen inherits a suture,
  basin, and Wilson-cycle parent.
- P82 also audits the current generated Earth-like Wilson objects and records
  the expected current residual: `tectonics.spreading_centers` are still absent
  even though ocean basins, Wilson cycles, gateways, rifts, passive margins,
  closing margins, and sutures exist.
- Current P82 output:
  `/Users/rayw/Projects/aevum/out_bench_p82_wilson_cycle_lifecycle_reference_20260627/`.
- Current P82 metrics:
  - scripted frame count: `7`;
  - scripted unique basin/lineage counts: `1/1`;
  - scripted object sets observed: `8/8`;
  - scripted gateway status count: `6`;
  - scripted parent link failures: `0`;
  - scripted old-orogen relicts: `1`;
  - scripted final basin age: `520 Myr`;
  - current ocean basin / Wilson cycle / gateway counts: `7/7/11`;
  - current rift/passive/spreading/closing/suture counts: `4/2/0/4/1`;
  - current basin stage count: `4`;
  - current active Wilson phase code count: `4`;
  - current missing object set: `tectonics.spreading_centers`.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p82_wilson_cycle_lifecycle_reference_pass -q`
  -> `1 passed in 31.14s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P82 --out out_bench_p82_wilson_cycle_lifecycle_reference_20260627`
  -> `status: pass`.

Next:

- Implement `P83.crust_sediment_province_coupling`.

2026-06-27 - Post-P82 staged research execution plan archived

- Updated `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md` with the detailed
  stage 0-7 execution control plan requested after the current real-Earth
  morphology review.
- The archived plan keeps this physiographic architecture track focused on
  real-Earth landform/process coverage: current gap inventory, source
  provenance, topography/bathymetry/planform envelopes, province-crust-
  sediment coupling, boundary/Wilson residual closure, mountain/rift/margin
  feature families, drainage/erosion/source-to-sink coupling, and integrated
  promotion audit.
- This document remains the province-architecture implementation log; the
  detailed source/theory/test/optimization checklist lives in the real-Earth
  research plan to avoid duplicating the same roadmap.

Next:

- Execute `P83.crust_sediment_province_coupling` as the first Stage 3
  implementation benchmark.

2026-06-27 - P83 crust-sediment-province coupling implemented

- Added `aevum/diagnostics/crust_sediment_province_coupling.py`.
- Added `P83.crust_sediment_province_coupling` and `run_p83_bench` to the
  tectonics benchmark CLI.
- P83 creates a deterministic synthetic continent fixture that couples
  province class, parent process, crust thickness, sediment thickness, basement
  age, stability, elevation, and relief.
- The fixture covers `13` province classes and `13` parent processes, and
  verifies that shields are old/stable without becoming default high flats,
  basins/passive margins are low by accommodation and sediment state, and
  orogens/LIPs are high because of parented crustal thickening.
- P83 also audits the current generated Earth-like world.  The current 900-cell
  audit has `31` continental landform objects, `8` landform kinds, low
  basin/lowland elevation relative to platforms, and high basin/lowland
  sediment relative to platforms.
- Recorded residual:
  production still lacks first-class `tectonics.continental_province_id`,
  `tectonics.continental_province_code`, and
  `tectonics.province_parent_process` fields.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p83_crust_sediment_province_coupling_20260627/`.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass -q`
  -> `1 passed in 30.85s`.
  Combined P79-P83 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p79_province_reference_graph_extraction_pass tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p81_boundary_process_geometry_reference_pass tests/test_tectonics_bench.py::test_p82_wilson_cycle_lifecycle_reference_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass -q`
  -> `5 passed in 123.70s`.

Next:

- Implement `P84.source_to_sink_sediment_budget`.

2026-06-27 - P84 source-to-sink sediment budget implemented

- Added `aevum/diagnostics/source_to_sink_sediment_budget.py`.
- Added `P84.source_to_sink_sediment_budget` and `run_p84_bench` to the
  tectonics benchmark CLI.
- P84 adds a deterministic sediment source-to-sink budget reference with
  mountain/platform sources and foreland, passive-margin, shelf, and ocean
  basin sinks.
- The fixture conserves source and sink sediment volume exactly at benchmark
  scale, keeps projected land/sea masks unchanged, keeps deposition below
  accommodation, and checks that routing edges close the zone budgets.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p84_source_to_sink_sediment_budget_20260627/`.
- Current P84 metrics:
  source/sink volumes `69000/69000 km3`, volume-balance fraction `0.0`, max
  accommodation utilization `0.68`, land-mask changes `0`, routing mismatches
  `0`.
- Current generated-world residual:
  `terrain.drainage_basins`, `terrain.sediment_routing_edges`, and
  `terrain.sediment_budget` are still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass -q`
  -> `1 passed in 31.08s`.
  P83-P84 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass -q`
  -> `2 passed in 61.70s`.

Next:

- Implement `P85.drainage_divide_province_alignment`.

2026-06-27 - P85 drainage-divide province alignment implemented

- Added `aevum/diagnostics/drainage_divide_province_alignment.py`.
- Added `P85.drainage_divide_province_alignment` and `run_p85_bench` to the
  tectonics benchmark CLI.
- P85 adds a deterministic multi-province drainage fixture with highland
  divides, west/interior, east/passive-margin, and south/rift drainage basins,
  plus explicit flow paths into plausible sinks.
- The reference gate checks divide/province-boundary alignment, highland
  alignment, flow-to-sink consistency, downhill paths, no divide/basin
  crossings, and contiguous non-checkerboard drainage basins.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p85_drainage_divide_province_alignment_20260627/`.
- Current P85 metrics:
  divide alignment `1.0`, highland alignment `1.0`, flow-to-sink consistency
  `1.0`, downhill-step fraction `1.0`, divide/basin crossings `0/0`, and max
  basin component count `1`.
- Current generated-world residual:
  `terrain.drainage_basins`, `terrain.drainage_divides`,
  `terrain.flow_direction`, and `terrain.flow_accumulation` are still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass -q`
  -> `1 passed in 31.11s`.
  P84-P85 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass -q`
  -> `2 passed in 61.83s`.

Next:

- Implement `P86.old_orogen_erosion_decay`.

2026-06-27 - P86 old-orogen erosion decay implemented

- Added `aevum/diagnostics/old_orogen_erosion_decay.py`.
- Added `P86.old_orogen_erosion_decay` and `run_p86_bench` to the tectonics
  benchmark CLI.
- P86 adds a deterministic old-orogen decay reference sequence from active
  collision orogen through post-collision high relief, decaying orogen,
  inherited old orogen, and subdued old-orogen province.
- The reference gate requires relief/elevation/crustal-root decay, late
  sediment-export decline after peak erosion, persistent inherited boundary
  trace and boundary strength, and parent-process links.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p86_old_orogen_erosion_decay_20260627/`.
- Current P86 metrics:
  relief `1900 -> 420 m`, relief-decay fraction `0.7789473684210526`,
  elevation `3200 -> 680 m`, boundary strength `1.0 -> 0.62`, minimum boundary
  overlap `1.0`, total sediment export `62160 km3`, final interval sediment
  export `9660 km3`, and parent-link failures `0`.
- Current generated-world audit:
  old-subdued-orogen objects `4`, parented old-subdued-orogen objects `4`,
  mean old-orogen elevation `912.553992412023 m`, mean old-orogen sediment
  thickness `881.3214900968968 m`, and old-orogen area
  `52719474237412.18 m2`.
- Current generated-world residual:
  `terrain.old_orogen_decay_stage`, `terrain.orogen_erosion_budget`,
  `terrain.orogen_boundary_memory`, and `terrain.orogen_sediment_export` are
  still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass -q`
  -> `1 passed in 30.75s`.
  P85-P86 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass -q`
  -> `2 passed in 61.75s`.

Next:

- Implement `P87.mountain_inventory_expression`.

2026-06-27 - P87 mountain inventory expression implemented

- Added `aevum/diagnostics/mountain_inventory_expression.py`.
- Added `P87.mountain_inventory_expression` and `run_p87_bench` to the
  tectonics benchmark CLI.
- P87 adds a GMBA-style small reference inventory for active margin ranges,
  collision ranges, broad collision plateaus, old subdued orogens,
  rift-shoulder ranges, volcanic arc chains, and extensional ranges.
- The reference gate requires object-backed ranges, valid hierarchy/parent
  range links, parent processes, finite area distribution, plausible elongation
  distribution, and mountain relief envelope.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p87_mountain_inventory_expression_20260627/`.
- Current P87 metrics:
  reference ranges/classes/levels `11/7/3`, reference total mountain area
  fraction `0.066`, reference median elongation `6.4`, and reference max relief
  `5000 m`.
- Current generated-world audit:
  mountain candidate objects `11`, expressed mountain objects `7`,
  parented/parent-linked mountain objects `11/11`, mountain kind count `2`,
  parent process/context coverage `1.0/1.0`, total mountain candidate area
  fraction `0.14999228221409389`, median elongation `1.8097984605150934`, and
  elongated mountain object count `0`.
- Current generated-world residual:
  `orogen` and `plateau` mountain kinds are still missing; first-class
  mountain inventory fields are still missing; elongated range expression is
  still underdeveloped.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass -q`
  -> `1 passed in 31.39s`.
  P86-P87 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass -q`
  -> `2 passed in 62.00s`.

Next:

- Implement `P88.rift_margin_escarpment_sequence`.

2026-06-27 - P88 rift-margin escarpment sequence implemented

- Added `aevum/diagnostics/rift_margin_escarpment_sequence.py`.
- Added `P88.rift_margin_escarpment_sequence` and `run_p88_bench` to the
  tectonics benchmark CLI.
- P88 adds a deterministic rift-to-passive-margin reference transect with
  platform, paired rift shoulders, rift basin/axis, passive-margin escarpment,
  coastal lowland, shelf, slope, rise, and abyssal plain.
- The reference gate requires adjacency and ordering across relief, sediment,
  shelf-slope-rise-abyss bathymetry, and parent-process context.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p88_rift_margin_escarpment_sequence_20260627/`.
- Current P88 metrics:
  reference zones/classes/edges `11/10/10`, missing reference edges `0`,
  escarpment relief `780 m`, shelf/slope/rise/abyss depths
  `120/1700/3100/4300 m`; current rift-basin/passive-lowland/wedge object
  counts `7/2/2`, lowland-near-shelf/wedge fractions `1.0/1.0`,
  rift-near-passive-margin fraction `0.985508126394939`, shelf depth p75
  `122.88691954044891 m`, and abyss depth p50 `3509.9321625103657 m`.
- Current generated-world residual:
  rift shoulder objects, escarpment objects, rift-margin sequence IDs/stages,
  and rift-margin lineage IDs are still missing; passive-margin lowland objects
  remain tiny.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass -q`
  -> `1 passed in 31.11s`.
  P87-P88 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass -q`
  -> `2 passed in 61.40s`.

Next:

- Implement `P89.plateau_area_cap_and_decay`.

2026-06-27 - P89 plateau area cap and decay implemented

- Added `aevum/diagnostics/plateau_area_cap_and_decay.py`.
- Added `P89.plateau_area_cap_and_decay` and `run_p89_bench` to the tectonics
  benchmark CLI.
- P89 gives the architecture track a reference for finite collision plateaus
  and volcanic/LIP plateaus with parent-process requirements, area caps, and
  decay.
- Current P89 metrics:
  reference variants/frames/stages `2/8/8`, max collision/volcanic plateau
  area fractions `0.024/0.016`, collision/volcanic elevation decay
  `2450/1050 m`, current plateau objects/detail cells `0/0`, current LIP
  objects `8`, and high interior without plateau support fraction
  `0.08958631689930298` of continental land.
- Current generated-world residual:
  first-class plateau inventory, age, decay stage, parent-process ID, lineage
  ID, plateau expression, and volcanic/LIP plateau expression are still
  missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p89_plateau_area_cap_and_decay_pass -q`
  -> `1 passed in 31.23s`.
  P88-P89 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass tests/test_tectonics_bench.py::test_p89_plateau_area_cap_and_decay_pass -q`
  -> `2 passed in 63.47s`.

Next:

- Implement `P90.current_world_morphology_gap_inventory`.

2026-06-27 - P90 current-world morphology gap inventory implemented

- Added `aevum/diagnostics/current_world_morphology_gap_inventory.py`.
- Added `P90.current_world_morphology_gap_inventory` and `run_p90_bench` to the
  tectonics benchmark CLI.
- P90 consolidates the current generated-world morphology residuals from
  P78-P89 plus compiler consistency into owner-layer and failure-category
  groups.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p90_current_world_morphology_gap_inventory_20260627/`.
- Current P90 metrics:
  gaps `58`, owner layers `8`, categories `8`, source suites `15`,
  unassigned/generic blockers `0/0`, current residual items `40`, missing
  review assets `8/8`, compiler passed envelope `True`, high-flat interior
  fraction `0.11796198797931971`, highland-without-parent fraction
  `0.12018884302629565`, basin/lowland fraction `0.391446377600517`, and
  major components `2`.
- Owner-layer counts:
  planform `9`, province graph `5`, boundary/lifecycle `2`, crust/sediment
  `1`, drainage/erosion `11`, landform expression `16`, bathymetry/margin `6`,
  and compiler/render `8`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/current_world_morphology_gap_inventory.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P90 --out out_bench_p90_current_world_morphology_gap_inventory_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p90_current_world_morphology_gap_inventory_pass -q`
  -> `1 passed in 32.57s`.
  P89-P90 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p89_plateau_area_cap_and_decay_pass tests/test_tectonics_bench.py::test_p90_current_world_morphology_gap_inventory_pass -q`
  -> `2 passed in 64.56s`.

Next:

- Implement `P91.integrated_real_earth_morphology_promotion_audit`.

2026-06-27 - P91 integrated real-Earth morphology promotion audit implemented

- Added `aevum/diagnostics/integrated_real_earth_morphology_promotion_audit.py`.
- Added `P91.integrated_real_earth_morphology_promotion_audit` and
  `run_p91_bench` to the tectonics benchmark CLI.
- P91 closes the stage 0-7 audit sequence: it checks P76-P90 stage summaries,
  archived 8000-cell P69 assets, fresh 900/2500-cell CI assets, compiler
  consistency, and a generated CI contact sheet.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p91_integrated_real_earth_morphology_promotion_audit_20260627/`.
- Current P91 metrics:
  stage pass `15/15`, high-resolution required PNGs `24/24`, CI required PNGs
  `16/16`, CI compilers `2/2`, CI inventories `2/2`, root P90 non-asset gaps
  `50`, owner layers `7`, and promotion blockers `9`.
- Promotion decision:
  audit completed, default Earth-like promotion remains off, and the next
  engineering entry is `P92.production_residual_owner_repair_plan`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_integrated_real_earth_morphology_promotion_audit_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p91_integrated_real_earth_morphology_promotion_audit_pass -q`
  -> `1 passed in 130.21s`.

Next:

- Define and implement `P92.production_residual_owner_repair_plan`.

2026-06-27 - P92 production residual owner repair plan implemented

- Added `aevum/diagnostics/production_residual_owner_repair_plan.py`.
- Added `P92.production_residual_owner_repair_plan` and `run_p92_bench` to the
  tectonics benchmark CLI.
- P92 assigns all P91 blockers and current residual items to ordered production
  repair packets, with explicit implementation targets, microbenchmarks,
  acceptance targets, validation suites, dependencies, and a final P91 reaudit.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p92_production_residual_owner_repair_plan_20260627/`.
- Current P92 metrics:
  blockers `9/9` assigned, owner layers `7/7` assigned, residual items
  `32/32` assigned, repair packets `8`, dependency order valid, climate/ocean/
  monsoon targets `0`, and final validation suite `P91`.
- Next implementation packet:
  `P92.1_planform_and_reference_calibration`; candidate suites are
  `P93.planform_reference_calibration` and
  `P93.generated_component_ribbon_envelope`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P92 --out out_bench_p92_production_residual_owner_repair_plan_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p92_production_residual_owner_repair_plan_pass -q`
  -> `1 passed in 0.87s`.

Next:

- Implement `P93.planform_reference_calibration`.

2026-06-27 - P93 planform reference calibration archived

- Added `aevum/diagnostics/planform_reference_calibration.py`.
- Added `P93.planform_reference_calibration`,
  `P93.generated_component_ribbon_envelope`, and `run_p93_bench` to the
  tectonics benchmark CLI.
- P93 consumes the archived P69/P78/P90/P91/P92 evidence and fixes the first
  production packet as a calibration problem, not as a visual cleanup pass.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p93_planform_reference_calibration_20260627/`.
- Current P93 metrics:
  P90 planform gaps `9`, calibration targets `5`, covered planform metrics
  `5`, cross-owner deferred targets `1`, P69 reference members `3`, P78
  current out-of-envelope count `5`, unresolved primary planform metrics `5`,
  and default promotion still blocked.
- Archived calibration directions:
  increase exposed land fraction, increase major land component count, reduce
  largest-component dominance where present, reduce exposed-land ribbon
  fraction, and reduce largest-landmass coastline over-complexity.  The
  trench-fraction excess is explicitly deferred to
  `P92.7_bathymetry_margin_sequence`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/planform_reference_calibration.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P93 --out out_bench_p93_planform_reference_calibration_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p93_planform_reference_calibration_pass -q`
  -> `1 passed in 0.79s`.

Next:

- Implement `P94.production_province_graph_fields` for
  `P92.2_production_province_graph_fields`.

2026-06-27 - P94 production province graph fields archived

- Added production continental province id/code/parent-process fields and
  object graphs under both terrain and tectonics namespaces.
- Added `aevum/diagnostics/production_province_graph.py` and the P94
  benchmark pair for field/object consistency and volcanic-LIP/rift adjacency.
- Updated P80/P83 expectations: production province graph is now available,
  and old province-graph pending residuals are cleared.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p94_production_province_graph_fields_20260627/`.
- Current P94 metrics:
  `2/2` microbenchmarks passed, both 900-cell and 2500-cell worlds ready,
  minimum province object count `43`, minimum class count `9`, minimum
  parent-process count `10`, maximum missing field/class/edge/object mismatch
  count `0`, and P80 900-cell missing reference classes/edges `0`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P94 --out out_bench_p94_production_province_graph_fields_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p94_production_province_graph_fields_pass -q`
  -> `3 passed in 179.99s`.

Next:

- Implement `P95.boundary_lifecycle_objects` for
  `P92.3_boundary_lifecycle_objects`.

2026-06-27 - Staged evidence collection archive cross-linked

- Cross-link entry retained for the P92/P93/P94/P95 sequence.
- Cross-linked the expanded evidence-collection archive in
  `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`.
- The continental physiographic architecture track should now use that archive
  as the source of truth for research batches, evidence-packet requirements,
  theory notes, benchmark enrichment, and optimization targets.

Next:

- Implement `P95.boundary_lifecycle_objects` for
  `P92.3_boundary_lifecycle_objects`.

2026-06-27 - P95 boundary lifecycle objects archived

- Added ridge-offset transform process geometry and sparse ridge/transform
  boundary-object aggregation.
- Updated P81/P82 expectations: transform process and spreading-center
  lifecycle objects are now available, and old residuals are cleared.
- Added the P95 benchmark pair for transform/spreading-center object coverage
  and current-world lifecycle audit.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p95_boundary_lifecycle_objects_20260627/`.
- Current P95 metrics:
  `2/2` microbenchmarks passed, both 900-cell and 2500-cell worlds ready,
  minimum transform cell count `15`, minimum ridge cell count `23`, minimum
  transform object count `3`, minimum ridge object count `4`, minimum
  spreading-center count `4`, missing process/object-set counts `0`, and
  parent link failures `0`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P95 --out out_bench_p95_boundary_lifecycle_objects_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p81_boundary_process_geometry_reference_pass tests/test_tectonics_bench.py::test_p82_wilson_cycle_lifecycle_reference_pass tests/test_tectonics_bench.py::test_p95_boundary_lifecycle_objects_pass -q`
  -> `3 passed in 177.16s`.

Next:

- Implement `P96.high_flat_interior_owner_reduction` and
  `P96.province_crust_sediment_surface_ordering` for
  `P92.4_crust_sediment_interior_relief_coupling`.

2026-06-27 - P96 crust/sediment interior relief coupling archived

- Added production surface ordering that couples continental province class,
  crustal support, sediment accommodation, passive-margin state, rift
  potential, and deformation context.
- Updated generated-world coupling audits to aggregate from
  `tectonics.continental_provinces` when the production graph is available,
  instead of treating legacy landform overlays as the current source of truth.
- Updated generated province-reference comparison so the production province
  graph is authoritative for class coverage and P80 residual clearance.
- Added/verified deterministic production coverage for platform, shield,
  passive-margin lowland, intracratonic basin, and foreland-basin anchors in
  generated major continents.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p96_crust_sediment_surface_ordering_20260627_d/`.
- Current P96 metrics:
  `2/2` microbenchmarks passed, both 900-cell and 2500-cell worlds ready,
  minimum production province class count `9`, minimum parent-process count
  `10`, maximum direct high-flat interior fraction `0.06246298636795882`,
  owner-aware P96 high-flat fraction after ordering `0.0`, minimum basin
  lowland fraction `0.37219087217151403`, minimum basin sediment excess over
  platforms `941.43280342115m`, and minimum orogen elevation excess over
  basins `554.1399306870508m`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/terrain.py aevum/diagnostics/crust_sediment_province_coupling.py aevum/diagnostics/generated_province_reference.py aevum/diagnostics/production_province_graph.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P96 --out out_bench_p96_crust_sediment_surface_ordering_20260627_d`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p94_production_province_graph_fields_pass tests/test_tectonics_bench.py::test_p96_crust_sediment_surface_ordering_pass -q`
  -> `4 passed in 306.66s`.

Next:

- Implement `P97.production_drainage_source_to_sink_fields` for
  `P92.5_drainage_source_to_sink_fields`.

2026-06-27 - P97 drainage/source-to-sink production fields archived

- Added production drainage, divide, flow, sediment routing, and old-orogen
  decay fields to terrain.
- Updated P84/P85/P86 generated-world audits so source-to-sink sediment budgets,
  drainage basin/divide alignment, and old-orogen erosion decay are validated
  against current production fields and objects.
- Added P97 generated-world diagnostics and connected P90/P91 reaudit logic to
  the current same-run stage summaries.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p97_drainage_source_to_sink_fields_20260627/`.
- Current P97 metrics:
  `2/2` microbenchmarks passed, both 900-cell and 2500-cell worlds ready,
  minimum drainage basin object count `18`, minimum routing edge count `39`,
  routing source/sink kind counts at least `3/5`, maximum sediment budget
  balance fraction `1.8436625312242433e-16`, maximum drainage divide fraction
  of land `0.5095258030376043`, divide alignment `1.0`, flow-to-sink
  consistency `1.0`, downhill path fraction `1.0`, and old-orogen decay budget
  fields present.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P97 --out out_bench_p97_drainage_source_to_sink_fields_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass tests/test_tectonics_bench.py::test_p97_drainage_source_to_sink_fields_pass -q`
  -> `4 passed in 228.06s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P90 --out out_bench_p90_after_p97_probe_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_after_p97_probe_20260627_b`
  -> `status: pass`.

Next:

- Implement `P98.production_landform_inventory_lifecycle_fields` for
  `P92.6_landform_inventory_lifecycle`.

2026-06-27 - P98 landform inventory/lifecycle production fields archived

- Added production mountain inventory fields and objects:
  `terrain.mountain_ranges`, `terrain.mountain_inventory`,
  `terrain.mountain_hierarchy_level`, `tectonics.mountain_belt_id`,
  `tectonics.mountain_parent_process_id`, and `terrain.mountain_ranges`
  objects.
- Added production plateau lifecycle fields and objects:
  `terrain.plateau_inventory`, `terrain.plateau_age_myr`,
  `terrain.plateau_decay_stage`, `terrain.plateau_parent_process_id`,
  `tectonics.plateau_lineage_id`, and `terrain.plateau_inventory` objects.
- Updated P87/P89 generated-world audits so mountain and plateau expression are
  validated from production fields, then updated P90/P91 same-run reaudit
  summaries.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p98_landform_inventory_lifecycle_20260627_d/`.
- Current P98 metrics:
  `2/2` microbenchmarks passed, both 900-cell and 2500-cell worlds ready,
  minimum production mountain range object count `6`, minimum mountain field id
  count `6`, minimum mountain inventory class count `4`, maximum missing
  mountain field/kind counts `0/0`, minimum plateau inventory cell count `19`,
  minimum volcanic/LIP plateau cells `19`, maximum missing plateau item/kind
  counts `0/0`, maximum total mountain area fraction `0.1244333061118279`, and
  maximum plateau area fraction `0.021075501103753724`.
- Reaudit:
  P90 after-P98 records `25` gaps, `13` current residual items, `4` owner
  layers, and no missing mountain/plateau production fields.
  P91 after-P98 records `17` non-asset root gaps, `3` root owner layers, `5`
  promotion blockers, and keeps default promotion blocked.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P98 --out out_bench_p98_landform_inventory_lifecycle_20260627_d`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P90 --out out_bench_p90_after_p98_probe_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_after_p98_probe_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py -k 'p87_mountain_inventory_expression or p89_plateau_area_cap_and_decay or p90_current_world_morphology_gap_inventory or p91_integrated_real_earth_morphology_promotion_audit or p98_landform_inventory_lifecycle'`
  -> `5 passed, 86 deselected, 2 warnings in 390.24s`.

Next at that point:

- Implement `P99.production_bathymetry_margin_sequence_fields` for
  `P92.7_bathymetry_margin_sequence`.

2026-06-27 - P99 bathymetry/margin sequence production fields archived

- Current output:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P99 --out out_bench_p99_bathymetry_margin_sequence_20260627_c`
  -> `status: pass`.
- Current P99 metrics:
  - `sequence_ready_world_count: 2`, `max_missing_sequence_item_count: 0`,
    `min_sequence_id_count: 1`, `min_lineage_id_count: 1`,
    `min_stage_count: 8`.
  - `min_rift_shoulder_cell_count: 143`,
    `min_escarpment_cell_count: 3`,
    `min_sequence_object_count: 1`.
  - `ordered_world_count: 2`, `min_shelf_stage_cell_count: 62`,
    `min_slope_stage_cell_count: 65`, `min_rise_stage_cell_count: 100`,
    `min_abyss_stage_cell_count: 119`.
- Current generated-world audit:
  - P88 now reports no missing rift-margin sequence items and no remaining
    rift-shoulder, escarpment, or rift-to-margin lineage residual.
  - P90 after-P99 records `19` gaps, `3` owner layers, `3` categories, and
    `8` current residual items, all of which are the required review assets.
  - P91 after-P99 records `11` non-asset root gaps, `2` root owner layers, and
    `4` named blockers: P69 reference calibration, P90 residuals,
    crust/sediment residuals, and planform residuals.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py -k 'p88_rift_margin_escarpment_sequence or p90_current_world_morphology_gap_inventory or p91_integrated_real_earth_morphology_promotion_audit or p99_bathymetry_margin_sequence'`
  -> `4 passed, 88 deselected, 2 warnings in 359.53s`.

Next at that point:

- Implement `P100.integrated_reaudit_and_promotion_gate` for
  `P92.8_integrated_reaudit_and_promotion_gate`.

2026-06-27 - P100 integrated reaudit and promotion gate archived

- Current output:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P100 --out out_bench_p100_integrated_reaudit_and_promotion_gate_20260627_c`
  -> `status: pass`.
- Current P100 metrics:
  - `repair_suite_count: 7`, `repair_suite_pass_count: 7`,
    `missing_repair_suite_count: 0`, `failing_repair_suite_count: 0`.
  - Fresh P91 reaudit: `p91_stage_suite_pass_count: 15`,
    `p91_ci_asset_set_complete_count: 2`,
    `p91_ci_compiler_passed_count: 2`.
  - Root residual state: `root_p90_non_asset_gap_count: 11`,
    `root_p90_owner_layer_count: 2`,
    `root_p90_residual_item_count: 0`.
  - Promotion blocker state: `promotion_blocker_count: 4`,
    `expected_blocker_set_matched: True`, `release_gate_allowed: False`.
- Current decision:
  - Default Earth-like promotion remains blocked.
  - P93-P99 repair suites have cleared the root blockers for province graph,
    boundary lifecycle, drainage/erosion, landform expression, and
    bathymetry/margin.
  - The remaining root blockers are planform and crust/sediment, plus the P69
    calibration and P90 residual umbrella blockers.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p100_integrated_reaudit_and_promotion_gate_pass -q`
  -> `1 passed, 2 warnings in 153.22s`.

Next:

- Define and implement `P101.planform_crust_sediment_residual_repair`.

2026-06-27 - P101+ real-Earth comparison repair plan archived

- Canonical plan:
  `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`, section
  `P101+ Real-Earth Comparison Repair Archive`.
- Current residual baseline:
  - `P100` passes with P93-P99 repair suites verified.
  - Default promotion remains blocked by P69 calibration, P90 current residuals,
    crust/sediment residuals, and planform residuals.
  - P90 after-P99 records `10` planform gaps, `1` crust/sediment gap, and `8`
    asset-review entries.
- Architecture implication:
  - The next work is not a broad new system and not a visual cleanup pass.
  - P101 should reproduce the residuals, attribute them to continent assembly,
    land exposure, sea-level solve, province/crust/sediment coupling, terrain
    base elevation, or compiler/render, and only then repair production behavior.
  - High flat interiors, low land fraction, overdominant largest component, low
    component count, excessive ribbon fraction, and coastline complexity should
    be treated as coupled planform plus crust/sediment failures.

Next:

- Start P101 with `P101.planform_residual_baseline` and
  `P101.crust_sediment_high_flat_repair` planning/diagnostics before changing
  terrain or tectonics production code.

2026-06-27 - P101 Phase 0 residual attribution archived

- Current output:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P101 --out out_bench_p101_phase0_current_residual_attribution_20260627`
  -> `status: pass`.
- Current architecture evidence:
  - Root 900-cell residuals are now reproduced and attributed: `10` planform
    gaps and `1` crust/sediment gap.
  - Planform defects are attributed to continent assembly, major component
    preservation, land exposure/sea-level solve, mainland ribbon pruning,
    coastline smoothing, and broad-interior support.
  - The high-flat interior defect is attributed to province/crust/sediment-driven
    interior elevation ordering, especially
    `terrain._regionalize_continental_surface`,
    `terrain._province_crust_sediment_surface_ordering`,
    `terrain._continental_detail_province`, and
    `terrain._production_continental_province_graph`.
  - 2500-cell cross-check confirms scale sensitivity: planform and
    crust/sediment remain present, and landform-expression residuals can reappear
    at higher resolution.  Future repairs must keep multi-resolution audits.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p101_current_residual_attribution_phase0_pass -q`
  -> `1 passed in 134.07s`.

Next:

- Phase 1 should enrich real-Earth reference/evidence packets and source
  fixtures for planform, province, crust, sediment, and multi-resolution
  comparisons.
- Phase 2 should start production planform repair only after the Phase 1 target
  envelopes and residual-owner risks are explicit.

2026-06-27 - P102 Phase 1 reference evidence packet matrix archived

- Current output:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P102 --out out_bench_p102_reference_evidence_packets_20260627`
  -> `status: pass`.
- Architecture evidence:
  - `R1_global_hypsometry_planform` binds planform/hypsometry/ocean metrics to
    P77/P78/P90/P101.
  - `R2_province_crust_sediment_basement` binds province, crust, sediment, and
    basement evidence to P79/P83/P96/P101.
  - `R3_boundary_wilson_deeptime` preserves boundary and deep-time lifecycle
    context for later planform repairs.
  - `R4_drainage_erosion_source_to_sink`, `R5_landform_margins_mountains_plateaus`,
    and `R6_case_study_feature_catalog` keep the cleared owner layers visible as
    regression risks.
  - The evidence packet matrix covers the P101 owner layers `planform` and
    `crust_sediment`, plus all already-cleared owner layers that must not regress.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p102_reference_evidence_packet_matrix_pass -q`
  -> `1 passed in 0.98s`.

Next:

- Start Phase 2/P103 planform mechanism repair.  The repair must improve land
  fraction, component count, largest-component dominance, ribbon fraction, and
  coastline complexity while preserving the P102 regression owners.

2026-06-27 - P103+ source-corpus enrichment plan archived

- Planning sync:
  - Canonical details live in
    `docs/REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md`, section
    `P103+ Source-Corpus Enrichment and Repair Planning Archive`.
  - The archive expands this physiographic architecture plan with staged source
    collection, real-Earth source tiers, theory-note requirements,
    microbenchmark reservations, and optimization discipline.
  - This is a planning-only archive update and does not change production
    terrain, tectonics, diagnostics, or thresholds.
- Execution implication:
  - P103 should repair planform mechanisms first, then later P104+ packets can
    address crust/sediment/interior elevation and natural province-boundary
    expression using the same evidence contract.
