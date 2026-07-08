# Real-Earth Geomorphology Research and Benchmark Plan

Status: active research archive
Owner: tectonics / terrain / diagnostics integration
Created: 2026-06-26
Scope: plate/continent/terrain work only; climate, ocean circulation, and
monsoon redesign stay paused except for regression fixes.

This document archives the staged research plan requested after the current
plate-engineering work reached the P70-P74 physiographic architecture layer.
The purpose is to turn real Earth topography and geomorphology into explicit
theory notes, metric targets, fixtures, implementation steps, and promotion
gates.  The target is not to copy modern Earth.  The target is that Aevum can
generate Earth-realistic landform classes for the right geological reasons,
with persistent objects and archive provenance rather than late visual repair.

Related documents:

- `docs/PLATE_TECTONICS_ENGINEERING_PLAN.md` - current P-series development
  log and release-gate history.
- `docs/PLATE_TECTONICS_REFACTOR_PLAN.md` - long-horizon replacement of
  random/heuristic major events with principle-model tectonics.
- `docs/CONTINENTAL_PHYSIOGRAPHIC_ARCHITECTURE_PLAN.md` - P70-P75 province
  architecture plan and implementation log.
- `docs/EARTH_GEOMORPHOLOGY_COVERAGE.md` - required Earth feature coverage
  contract.

## Current Engineering State

As of 2026-06-27:

- P70 has an offline research-source and metric-schema scaffold.
- P71 has deterministic physiographic province graph fixtures.
- P72 has a generated-world province diversity gate.
- P73 has offline real-Earth case-study sketches for major continents.
- P74 has production terrain coupling for province templates, including
  `terrain.passive_margin_lowland` and passive-margin lowland landform objects.
- P75 has completed the release/promotion audit.  It verifies P69 assets,
  P70-P74 gates, and legacy P29/P48/P49/P68 gates, then records that default
  promotion must remain blocked until Stage A reference extraction and
  first-class production province graph work are done.
- P76-P78 have converted source-ledger and real-Earth hypsometry/bathymetry
  targets into executable offline diagnostics.
- P79-P80 now compare generated major-continent province graphs against the
  derived real-Earth province reference graph.
- P81-P82 now provide executable boundary-process and Wilson-cycle lifecycle
  references, with current generated-world residuals recorded for transform
  boundaries and spreading-center objects.
- P83-P89 now provide executable fixtures and generated-world audits for
  crust/sediment/province coupling, source-to-sink sediment budgets, drainage
  divides, old-orogen erosion decay, mountain inventories, and
  rift-to-passive-margin sequences, plus plateau area caps and lifecycle decay.
- P90 now provides a current-world morphology gap inventory.  It groups the
  generated world's remaining defects by owner layer and failure category, and
  confirms that the next implementation entry is an integrated P91 promotion
  audit rather than another ad hoc visual pass.
- P91 now provides the integrated real-Earth morphology promotion audit.  It
  verifies P76-P90 stage summaries, archived 8000-cell P69 PNG evidence, fresh
  900/2500-cell CI world assets, compiler consistency, and contact-sheet
  availability.  The audit is complete but keeps default Earth-like promotion
  blocked behind named residual owner layers.
- P92 now provides the production residual owner repair plan.  It assigns all
  P91 promotion blockers and all current residual items to ordered production
  repair packets with implementation targets, microbenchmarks, acceptance
  targets, validation suites, dependencies, and a final P91 reaudit gate.
- P93 now archives the planform reference calibration entry.  It binds the
  P69/P78/P90/P91/P92 evidence into explicit component, ribbon, coastline, land
  fraction, and largest-component calibration targets while keeping default
  promotion blocked and deferring trench over-expression to the later
  bathymetry/margin packet.
- P94 now provides the production continental province graph fields and
  objects, including required province classes and reference graph edges.
- P95 now provides production boundary lifecycle objects for transforms,
  ridges, spreading centers, and Wilson-cycle current-world readiness.
- P96 now provides production crust/sediment/interior-relief coupling.  Current
  coupling audits aggregate from `tectonics.continental_provinces`, generated
  province-reference comparison treats the production province graph as
  authoritative, and surface ordering now links basin/platform/orogen elevation
  and sediment signals to crust and province state.
- P97 now provides production drainage/source-to-sink fields and objects,
  including drainage basins, divides, flow direction/accumulation, sediment
  source/sink budgets, routing edges, and old-orogen erosion-decay fields.
- P98 now provides production landform inventory/lifecycle fields and objects
  for mountain ranges and plateaus, including hierarchy/parent-process fields,
  plateau age/decay stage, LIP/volcanic plateau expression, and P90/P91
  after-P98 reaudit evidence.
- P99 now provides production bathymetry/margin sequence fields and objects,
  including rift shoulders, escarpments, rift-to-passive-margin sequence
  ids/stages, lineage ids, and auditable shelf/slope/rise/abyss ordering.
- P100 now provides the integrated after-P99 reaudit and promotion decision
  gate.  It verifies P93-P99 repair-suite evidence, reruns P91, confirms
  default promotion must remain closed, and narrows the next production work to
  planform plus crust/sediment residual repair.
- P101 Phase 0 now reproduces and attributes the remaining planform plus
  crust/sediment residuals without changing production behavior.
- P102 now archives reference evidence packets that bind sources, theory
  claims, derived metrics, generated-world audits, optimization targets, and
  residual owner layers for the next repairs.
- The P103+ source-corpus enrichment archive now defines how to collect
  substantial real-Earth evidence and convert it into benchmark details,
  theory notes, implementation methods, and optimization targets before later
  production tuning.

Known residual blocker:

- Generated continents can now expose object-backed province, boundary,
  crust/sediment, drainage/source-to-sink, mountain, plateau, and
  bathymetry/margin sequence structure.  The remaining blockers are now
  concentrated in planform and crust/sediment calibration, plus the P69
  Earth-like reference calibration and P90 residual umbrella blockers.
- The next implementation entries are:
  `P103.planform_mechanism_repair` for remaining land/planform residuals, then
  `P104A.continental_mosaic_object_expression` before any broad P104 elevation
  retuning.  The current map review shows that large continents can still look
  internally monotonous because production province objects are richer than the
  rendered/elevation response that consumes them.

## Research Principles

- Use real Earth as a feature-class and distribution reference, not as a map to
  clone.
- Treat source data as evidence for process relationships: province adjacency,
  relief envelopes, crustal setting, sediment routing, drainage geometry, and
  lifecycle persistence.
- Keep large data out of the repo.  Store source ledgers, extraction scripts,
  and small derived metric fixtures.
- Benchmarks must run offline.  Data acquisition and raster/vector extraction
  are separate preparation steps.
- A feature is covered only when it has a parent process object, a lifecycle
  state, generated-world expression, and metric evidence.
- Random seeds may perturb textures and ensembles; they must not choose major
  geological events, province existence, or final seaway openings.

## Post-P82 Execution Control Plan

This section is the archived plan for the next research and implementation
cycle.  It responds to the current priority: compare generated landforms with
real Earth first, repair missing geomorphic/tectonic mechanisms, and keep
climate, ocean circulation, and monsoon work paused until the land/terrain
architecture is credible.

Execution rule for every stage:

- collect or reference sources before changing production behavior;
- write the theory note as process relationships, not visual descriptions;
- create a deterministic fixture or generated-world diagnostic before tuning;
- record residual blockers explicitly instead of hiding them in permissive
  gates;
- update this document and the P-series engineering log after each benchmark.

### Stage 0. Current-State Gap Inventory

Purpose:

- Freeze what the current generator can and cannot express after P89.
- Convert visual complaints into named gap categories: missing process, wrong
  amplitude, wrong area scale, wrong adjacency, wrong lifecycle, wrong
  sediment/crust coupling, and compiler/render mismatch.

Inputs:

- P69 high-resolution audit assets;
- P78 generated hypsometry comparison;
- P80 generated/reference province graph comparison;
- P81 boundary-process geometry audit;
- P82 Wilson-cycle lifecycle audit;
- current `elevation.png`, `terrain_provinces.png`,
  `continental_detail_provinces.png`, `ocean_depth_provinces.png`, and
  `hexmap.png` asset sets.

Theory questions:

- Which defects are failures of continent planform, and which are failures of
  province expression on otherwise acceptable continents?
- Which defects can be fixed by adding parent process objects, and which need
  archive/lifecycle continuity?
- Which apparent terrain problems are really missing crust, sediment, drainage,
  or erosion state?

Implementation method:

- Add a compact gap-inventory diagnostic only if current P78/P80/P81/P82
  metrics cannot identify the failing layer.
- Keep this stage read-only for production code.

Microbenchmark target:

- Candidate `P90.current_world_morphology_gap_inventory`.
- Metrics: current residual categories, missing object fields, missing map
  assets, high-flat interior share, unparented highland share, basin/lowland
  share, and rendered/compiler mismatch rates.
- Acceptance: every residual maps to a named future stage; no generic
  "looks wrong" blocker remains.

### Stage 1. Source Corpus and Extraction Provenance

Purpose:

- Expand P76 from metadata coverage into reproducible source-acquisition and
  extraction instructions, while keeping raw rasters/vectors outside git.

Source groups:

- ETOPO/GEBCO/Natural Earth for topography, bathymetry, coastlines, shelves,
  islands, and planform.
- USGS/NPS/global tectonic province sources for province boundaries and
  reference class graphs.
- PB2002/GPlates/EarthByte/GEM/World Stress Map for plate boundaries, faults,
  stress regimes, and deep-time case studies.
- CRUST1.0/GLiM/NOAA sediment/GlobSed for crustal thickness, lithology, and
  sediment thickness.
- HydroSHEDS/HydroBASINS/HydroRIVERS and GMBA for drainage, erosion context,
  and mountain inventories.

Theory questions:

- Which sources define first-order process relationships, and which only
  provide calibration envelopes?
- What can be distributed as small derived fixtures, and what must remain a
  local regeneration input?

Implementation method:

- Add extractor entry points that write small JSON metric fixtures with source
  IDs, version notes, extraction parameters, projection notes, and checksums.
- Do not make runtime generation depend on web access or raw external data.

Microbenchmark target:

- Extend `P76.reference_source_ledger_schema` rather than replacing it.
- Metrics: extraction scripts registered, derived fixture checksums present,
  raw-data storage policy explicit, and source IDs stable across fixtures.
- Acceptance: a new developer can regenerate every small fixture from the
  ledger without guessing source versions or projection choices.

### Stage 2. Global Hypsometry, Bathymetry, and Planform Envelopes

Purpose:

- Tighten P77/P78 so global elevation, bathymetry, shelves, trenches, islands,
  ribbons, and coastline complexity are calibrated against real Earth feature
  classes rather than broad visual expectations.

Theory questions:

- What is the acceptable envelope for Earth-like land fraction, major component
  count, largest component share, island area distribution, and narrow-neck
  density?
- How much high terrain should exist globally and per major continent?
- How should continental shelves, slopes, rises, abyssal plains, ridges, and
  trenches distribute with distance from coast and plate boundaries?

Implementation method:

- Add equal-area metric extraction for planform and coastline classes.
- Compare generated worlds at CI scale first, then confirm at 8000-cell audit
  scale before changing promotion status.

Microbenchmark targets:

- Extend `P77.real_earth_hypsometry_extraction` with planform bins.
- Extend `P78.generated_hypsometry_envelope` with component, island, ribbon,
  coastline, shelf, ridge, trench, and abyssal-fraction diagnostics.
- Acceptance: basic hypsometry stays green while failures identify whether the
  cause is planform generation, bathymetry generation, or map compilation.

### Stage 3. Province, Crust, Sediment, and Basement Coupling

Purpose:

- Resolve the core inland problem: continents should contain multiple natural
  physiographic provinces, and interior elevation should follow crustal
  support, sediment accommodation, basement age, lithology, and inherited
  structure rather than becoming uniformly high and flat.

Theory questions:

- Why are shields old and stable without being default high plateaus?
- How do covered platforms, intracratonic basins, foreland basins, failed
  rifts, and passive-margin lowlands become low without erasing continental
  identity?
- How should crustal thickness and sediment thickness interact with base
  elevation and local relief?

Implementation method:

- Add first-class production province IDs/codes and parent links where the
  current raster/detail fields are only proxies.
- Let terrain base elevation read province template, crust thickness,
  basement age, sediment accommodation, and parent process before final
  texture.
- Keep sediment routing conservative until source-to-sink budget tests pass.

Microbenchmark targets:

- `P83.crust_sediment_province_coupling`.
- `P84.source_to_sink_sediment_budget`.
- Follow-up generated-world gate for per-continent province/crust/sediment
  ordering after the synthetic fixture passes.
- Acceptance: basins/passive margins are low by accommodation and sediment
  state; shields remain old/stable but not uniform high flats; major
  highlands/plateaus have parent objects and area limits.

### Stage 4. Boundary Geometry and Deep-Time Lifecycle Completion

Purpose:

- Close the residuals recorded by P81/P82 before using boundary objects as
  stronger terrain parents.

Theory questions:

- How should transform offsets, ridge segmentation, fracture zones, trenches,
  sutures, diffuse deformation belts, and passive margins coexist on a sphere?
- How should a Wilson-cycle archive preserve basin, margin, gateway, suture,
  and old-orogen identity after active motion fades?

Implementation method:

- Add first-class transform/spreading-center object expression from existing
  boundary masks and ocean-basin archive state.
- Preserve spherical continuity across the antimeridian and avoid rectangular
  render seams in geometry logic.

Microbenchmark targets:

- Successors to `P81.boundary_process_geometry_reference` and
  `P82.wilson_cycle_lifecycle_reference`.
- Metrics: transform boundary cells/objects, spreading-center objects,
  ridge-transform adjacency, persistent basin lineage, and inherited suture to
  old-orogen links.
- Acceptance: current residuals `transform` and
  `tectonics.spreading_centers` are removed or replaced by narrower recorded
  blockers.

### Stage 5. Mountains, Plateaus, Rifts, Margins, and Special Landforms

Purpose:

- Ensure Earth's major landform families can be generated by process models:
  active orogens, old orogens, collision plateaus, volcanic/LIP plateaus,
  rift systems, passive-margin escarpments, arcs, forearcs, back-arcs,
  hotspot chains, marginal seas, and ocean gateways.

Theory questions:

- Which mountain belts should be narrow ranges, broad plateaus, old subdued
  belts, or arc chains?
- How large can collision and volcanic plateaus be before they become the
  default interior state?
- How do rift basins, shoulders, escarpments, coastal plains, shelves, and
  sediment wedges form an ordered margin sequence?

Implementation method:

- Promote feature objects from classification results into terrain drivers.
- Add area caps, lifecycle decay, and parent-process requirements for each
  high-relief or special landform class.

Microbenchmark targets:

- `P87.mountain_inventory_expression`.
- `P88.rift_margin_escarpment_sequence`.
- `P89.plateau_area_cap_and_decay`.
- Acceptance: mountain and plateau features are object-backed, finite in area,
  and correctly adjacent to basins, margins, rifts, arcs, or sutures.

### Stage 6. Drainage, Erosion, and Source-to-Sink Surface Processes

Purpose:

- Move final terrain expression from cellwise response toward region, basin,
  and catchment response.

Theory questions:

- How should drainage divides follow highlands, shields, sutures, rift
  shoulders, and passive-margin escarpments?
- How does erosion lower old orogens and mature interiors while preserving
  boundaries?
- How should sediment move from highland sources to foreland, intracratonic,
  passive-margin, delta/fan, shelf, and abyssal sinks?

Implementation method:

- Build drainage and sediment paths from province/base-relief surfaces before
  small-scale terrain texture.
- Keep climate-dependent discharge as a placeholder until climate redesign
  resumes.

Microbenchmark targets:

- `P85.drainage_divide_province_alignment`.
- `P86.old_orogen_erosion_decay`.
- Extended `P84.source_to_sink_sediment_budget`.
- Acceptance: drainage paths are coherent at region scale; sediment volume is
  conserved within tolerance; erosion changes relief without breaking land/sea
  masks.

### Stage 7. Integrated Earth-Like Promotion Audit

Purpose:

- Decide whether the new terrain/plate architecture is good enough to replace
  the current flagged Earth-like behavior.

Inputs:

- Updated P76-P89 diagnostics;
- 900/2500-cell CI worlds;
- 8000-cell high-resolution reference assets;
- optional higher-resolution profile only after performance evidence.

Implementation method:

- Generate contact sheets for elevation, province IDs/classes, crust age,
  crust thickness, sediment thickness, drainage, landforms, ocean provinces,
  boundary processes, history/timeline, and hexmap.
- Promote only after metrics and visual assets agree.

Microbenchmark target:

- Candidate `P91.integrated_real_earth_morphology_promotion_audit`.
- Metrics: pass counts across P76-P90, current residual list, key PNG asset
  availability, component/planform metrics, province/crust/sediment ordering,
  drainage/sediment budget, and compiler consistency.
- Acceptance: generated Earth-like worlds can express real-Earth feature
  classes for geological reasons, while wet/dry variants remain plausible
  alternatives rather than forced modern-Earth copies.

## Source Corpus To Build

### S0. Source Ledger and Reproducibility

Goal:

- Create a structured source ledger that records URL, version, license/use
  notes, downloaded file names, extraction date, projection, resolution, and
  derived-metric files.

First implementation:

- Add a small JSON/YAML source ledger under `data/reference/` only after the
  first extraction script exists.
- Keep raw downloaded rasters/vectors outside git by default.
- Every derived metric fixture must include `source_id`, source version,
  extraction code hash, and preprocessing steps.

Microbenchmarks:

- `P76.reference_source_ledger_schema`
  - fixture: local source ledger with no raw data requirement;
  - metrics: required metadata fields present, source IDs stable, license
    status explicit;
  - acceptance: no source without version/license/projection notes.

### S1. Global Topography, Bathymetry, and Planform

Primary sources:

- NOAA ETOPO Global Relief Model:
  https://www.ncei.noaa.gov/products/etopo-global-relief-model
- NOAA ETOPO 2022 metadata:
  https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem%3Aetopo_2022
- GEBCO Gridded Bathymetry Data:
  https://www.gebco.net/data-products/gridded-bathymetry-data
- Natural Earth 1:10m physical vectors:
  https://www.naturalearthdata.com/downloads/10m-physical-vectors/

Questions:

- What are the global land and ocean hypsometric envelopes?
- How much land sits in lowlands, moderate uplands, high plateaus, and high
  mountain terrain?
- How wide are shelves, slopes, abyssal plains, trenches, and mid-ocean-ridge
  highs relative to coastlines and plate boundaries?
- What are plausible component counts, coastline complexity ranges, island
  fractions, and narrow-neck/ribbon fractions at Aevum's working resolutions?

Derived metrics:

- land/ocean hypsometry quantiles and histogram bins;
- land fraction and major land component share;
- lowland fractions below 500 m and 1000 m;
- highland fractions above 1500 m, 2500 m, and 4000 m;
- local relief quantiles by moving-window radius;
- coastline complexity by component size;
- island-count distribution by area bin;
- shelf/slope/rise/abyss/trench fractions and nearshore depth envelope;
- ridge-high and transform/fracture-zone lineament density.

Implementation notes:

- Start with coarse downsampled global rasters so tests stay fast.
- Use equal-area or cell-area weighted statistics; do not trust raw equirect
  pixel counts near the poles.
- Store only small extracted metric JSON fixtures, not full ETOPO/GEBCO rasters.

Microbenchmarks:

- `P77.real_earth_hypsometry_extraction`
  - fixture: downsampled ETOPO/GEBCO metric file;
  - metrics: land/ocean hypsometry, shelf/slope/abyss fractions, highland tail;
  - acceptance: extraction reproduces stored checksum and expected broad
    values within tolerance.
- `P78.generated_hypsometry_envelope`
  - fixture: Earth-like generated worlds at 2500/8000 cells;
  - metrics: same schema as `P77`;
  - acceptance: generated distributions are inside initial Earth-like envelope
    or explicitly marked as designed variant behavior.

### S2. Physiographic and Geological Provinces

Primary sources:

- USGS Physiographic Divisions of the conterminous United States:
  https://data.usgs.gov/datacatalog/data/USGS%3Ae04ea9e9-17b6-45ae-b279-7bc35ea79539
- NPS physiographic province explanation:
  https://www.nps.gov/subjects/geology/physiographic-provinces.htm
- Global geologic provinces and tectonic plates:
  https://zenodo.org/records/6586972
- Hasterok global tectonics repository:
  https://github.com/dhasterok/global_tectonics

Questions:

- Which province classes are needed for generated worlds to read as natural
  continents?
- Which province boundaries should persist after active tectonics fades?
- What adjacency graphs are common: shield-platform-basin, orogen-foreland,
  rift-shoulder-basin, passive-margin-lowland-shelf, arc-forearc-trench?
- How many distinct provinces should a large continent contain at Aevum
  resolution?

Derived metrics:

- province class shares by continent and by case-study region;
- province count per major component;
- largest province fraction per major component;
- adjacency edge inventory and forbidden adjacency list;
- boundary length and boundary-strength distribution;
- province relief envelope: base elevation, local relief, slope, roughness;
- parent process coverage for each province.

Implementation notes:

- Keep P71 fixtures as the first synthetic grammar.
- Build reference sketches first; only then add raster/vector extraction.
- Do not require generated continents to match North America or Africa.  Require
  the same feature-class graph to be generatable.

Microbenchmarks:

- `P79.province_reference_graph_extraction`
  - fixture: extracted USGS/NPS/global-province reference graph metrics;
  - metrics: class shares, adjacency edges, largest province fraction;
  - acceptance: extracted graph is deterministic and complete enough for
    calibration.
- `P80.generated_major_continent_province_graph`
  - fixture: Earth-like generated worlds;
  - metrics: per-continent class count, adjacency, parent process coverage;
  - acceptance: every large continent has multiple process-parented provinces,
    no single smooth platform dominates unless the continent is small.

### S3. Plate Kinematics, Boundaries, and Deep-Time Context

Primary sources:

- GPlates:
  https://www.gplates.org/
- EarthByte data and models:
  https://www.earthbyte.org/category/resources/data-models/
- GPlates Portal:
  https://portal.gplates.org/
- PB2002 plate boundary model:
  https://peterbird.name/publications/2003_pb2002/2003_pb2002.htm
- World Stress Map:
  https://www.world-stress-map.org/
- GEM Global Active Faults:
  https://github.com/GEMScienceTools/gem-global-active-faults

Questions:

- Which present-day boundary geometries and stress regimes should be reflected
  by object rules?
- How should deep-time reconstructions inform Wilson-cycle object lifecycles
  without forcing Aevum to replay named Earth history?
- How should inherited sutures, old orogens, rifts, and failed arms constrain
  later deformation?

Derived metrics:

- boundary-type length fractions: ridge, transform, subduction, collision,
  diffuse deformation;
- age/state transitions for rift -> passive margin -> ocean basin -> active
  margin -> suture;
- active-fault density in deformation belts versus stable interiors;
- stress-regime consistency for generated rift/collision/shear zones.

Implementation notes:

- Use GPlates/EarthByte first as case-study evidence and future tooling, not as
  a hard runtime dependency.
- Present-day boundary maps should calibrate geometry and relationships; they
  should not become a literal template for every generated world.

Microbenchmarks:

- `P81.boundary_process_geometry_reference`
  - fixture: present-day plate boundary reference and synthetic generated
    boundary networks;
  - metrics: boundary length shares, curvature, transform offsets, trench/ridge
    adjacency;
  - acceptance: generated boundary networks include plausible type diversity and
    spherical continuity.
- `P82.wilson_cycle_lifecycle_reference`
  - fixture: scripted cycle with rift, ocean opening, passive margin,
    subduction, collision, suture, old orogen;
  - metrics: lifecycle object survival, archive continuity, province
    inheritance;
  - acceptance: no stage depends only on current-frame random labels.

### S4. Crust, Lithology, Sediment, and Basement State

Primary sources:

- CRUST1.0:
  https://ds.iris.edu/ds/products/emc-crust10/
- GLiM Global Lithological Map:
  https://www.geo.uni-hamburg.de/en/geologie/forschung/aquatische-geochemie/glim.html
- GLiM paper/data overview:
  https://doi.pangaea.de/10.1594/PANGAEA.788537
- NOAA total sediment thickness:
  https://www.ncei.noaa.gov/products/total-sediment-thickness-oceans-seas
- GlobSed metadata:
  https://data.noaa.gov/onestop/collections/details/0ed0104f-2add-4c87-8fb3-9787b6d416c7

Questions:

- How should crustal thickness, basement age, sediment thickness, and lithology
  jointly determine platform, shield, basin, orogen, and margin expression?
- Why should stable interiors not become uniformly high and flat?
- How much sediment should accumulate in intracratonic, foreland, passive
  margin, delta/fan, and abyssal contexts?

Derived metrics:

- crustal thickness distributions by province class;
- sediment thickness by basin/margin/ocean setting;
- shield/platform/basin lithology class mix;
- basement-age and rework-state proxies by province;
- expected base-elevation envelope by crustal thickness and province support.

Implementation notes:

- Production terrain should read an explicit basement/province state before it
  uses final elevation thresholds.
- Sediment routing must conserve source-to-sink volume at benchmark scale
  before it is allowed to alter production elevation.

Microbenchmarks:

- `P83.crust_sediment_province_coupling`
  - fixture: shield/platform/basin/orogen/passive-margin synthetic continent;
  - metrics: crust thickness, sediment, elevation, relief ordering;
  - acceptance: basins and passive margins are low without erasing parent
    continent; shields are old/stable but not necessarily high and flat.
- `P84.source_to_sink_sediment_budget`
  - fixture: mountain source, foreland sink, passive-margin sink, ocean basin;
  - metrics: sediment mass balance and accommodation;
  - acceptance: no land-mask regression and conserved sediment within tolerance.

### S5. Drainage, Erosion, and Surface Process Coupling

Primary sources:

- HydroSHEDS:
  https://www.hydrosheds.org/
- HydroBASINS:
  https://www.hydrosheds.org/products/hydrobasins
- HydroRIVERS:
  https://www.hydrosheds.org/products/hydrorivers

Questions:

- How should drainage divides follow mountain belts, plateaus, shields, rift
  shoulders, old sutures, and passive-margin escarpments?
- How should erosion lower old orogens and mature platforms without destroying
  active orogens and plateaus too quickly?
- How should drainage basins route sediment toward foreland, intracratonic,
  passive-margin, delta/fan, and abyssal sinks?

Derived metrics:

- drainage basin count and size distribution;
- main-river length and branching hierarchy;
- basin hypsometric integral;
- divide alignment with province boundaries and highlands;
- sediment export potential by source province;
- delta/fan candidate count and margin adjacency.

Implementation notes:

- Drainage should be generated from province/base-relief surfaces before final
  texture and small noise.
- Climate-dependent discharge can be a placeholder for now, because climate
  redesign is paused; use topography-only drainage structure and conservative
  sediment capacity.

Microbenchmarks:

- `P85.drainage_divide_province_alignment`
  - fixture: multi-province continent with orogen, shield, basin, passive
    margin, rift shoulder;
  - metrics: divide-boundary alignment, flow-to-sink consistency;
  - acceptance: basins drain to plausible sinks and divides avoid random
    checkerboard paths.
- `P86.old_orogen_erosion_decay`
  - fixture: active orogen aging into old orogen and subdued belt;
  - metrics: relief decay, sediment output, boundary persistence;
  - acceptance: old orogen remains a province boundary after relief decays.

### S6. Mountains, Plateaus, Rifts, Margins, and Special Landforms

Primary sources:

- GMBA Mountain Inventory:
  https://www.earthenv.org/mountains
- GMBA mountain inventory paper:
  https://www.nature.com/articles/s41597-022-01256-y
- Global geologic provinces and tectonic plates:
  https://zenodo.org/records/6586972

Case-study families:

- active collision belts and broad plateaus;
- active continental margins and volcanic arcs;
- old subdued orogens and suture belts;
- intracontinental rifts and failed rifts;
- passive-margin coastal plains and escarpments;
- LIP/volcanic plateau surfaces;
- hotspot island/seamount chains;
- marginal seas and ocean gateways.

Derived metrics:

- mountain range count and hierarchy by continent;
- range area distribution and elongation;
- relief amplitude by mountain type;
- plateau area caps and edge steepness;
- rift-basin length/width and shoulder asymmetry;
- passive-margin lowland width and shelf adjacency;
- hotspot/LIP age progression and finite area.

Microbenchmarks:

- `P87.mountain_inventory_expression`
  - fixture: GMBA-style reference metrics plus generated mountain objects;
  - metrics: range count, hierarchy, elongation, parent process coverage;
  - acceptance: mountains are object-backed and not merely thresholded noise.
- `P88.rift_margin_escarpment_sequence`
  - fixture: rift opening to passive margin with escarpment and shelf;
  - metrics: rift basin, shoulder, escarpment, lowland, shelf adjacency;
  - acceptance: margin lowland and shelf are coupled without deep nearshore
    pits or far-ocean shoals.
- `P89.plateau_area_cap_and_decay`
  - fixture: collision plateau and volcanic plateau variants;
  - metrics: plateau area, relief, parent object, decay;
  - acceptance: plateaus are present but area-limited and not default interiors.

## Staged Evidence Collection and Benchmark Enrichment Archive

Archived: 2026-06-27

This section turns the research request into an execution archive.  The goal is
to collect enough real-Earth evidence to enrich four things before further
large production rewrites:

- test details: what a benchmark must measure and what failure means;
- theory basis: which process relationship the metric is defending;
- implementation method: which object, field, or archive layer should own the
  behavior;
- optimization target: which parameter or rule can be tuned, and which result
  should be considered invalid rather than merely ugly.

### Evidence Packet Contract

Every collected source batch should produce a small evidence packet.  Raw
rasters/vectors stay outside git unless they are explicitly tiny fixtures.
The repo should keep only extraction code, source metadata, and derived
metrics.

Required packet fields:

- `packet_id`: stable id such as `R1_global_hypsometry_planform`.
- `source_ids`: source names and versions from the P76 ledger.
- `raw_data_policy`: local path expectation, license note, and no-raw-data
  repo policy.
- `extraction_method`: projection, downsampling, masking, equal-area weighting,
  and checksum strategy.
- `theory_claims`: process relationships the packet supports.
- `derived_metrics`: JSON-compatible metric names and units.
- `reference_fixture`: synthetic or derived small fixture used by tests.
- `generated_world_audit`: how current Earth-like worlds are measured with the
  same schema.
- `optimization_targets`: envelope values or monotonic relationships.
- `residual_policy`: which failures block promotion now and which remain
  recorded future work.
- `asset_review`: PNG/contact-sheet outputs required for human inspection.

Microbenchmarks created from these packets must have two tracks:

- a reference track that proves the extracted or hand-curated evidence is
  internally coherent;
- a generated-world track that audits current output and records residuals
  without weakening the reference standard.

### R0. Current Generated-World Forensic Baseline

Purpose:

- Freeze the current artifacts before more tuning.
- Explain current failures by layer: planform, province graph, crust/sediment
  support, drainage/erosion, landform object coverage, bathymetry, or compiler.

Inputs:

- P69 high-resolution physical ensemble assets.
- P76-P89 benchmark outputs.
- Current Earth-like `elevation.png`, `terrain_provinces.png`,
  `continental_detail_provinces.png`, `ocean_depth_provinces.png`,
  `crust_age.png`, `history.png`, `timeline.png`, and `hexmap.png`.

Metrics to add or tighten:

- high-flat interior share by continent and by province;
- highland cells lacking active/old-orogen, plateau, rift-shoulder, arc, or
  plume parentage;
- per-continent lowland share below 500 m and 1000 m;
- per-continent province count and largest province fraction;
- shelf/nearshore depth envelope and far-ocean shoal fraction;
- compiler/elevation sign mismatch and class mismatch;
- current residual list grouped by missing object, missing field, wrong
  amplitude, wrong adjacency, wrong lifecycle, or rendering mismatch.

Candidate benchmark:

- `P90.current_world_morphology_gap_inventory`.

Exit criteria:

- Every visible defect has a named owner stage and no generic visual complaint
  remains as an implementation instruction.

### R1. Global Topography, Bathymetry, and Planform Evidence

Purpose:

- Anchor global land/ocean geometry, elevation, and depth envelopes to real
  Earth feature classes.
- Separate acceptable Earth-like variation from invalid morphology.

Sources:

- NOAA ETOPO 2022 Global Relief Model.
- GEBCO Gridded Bathymetry.
- Natural Earth land/coast/island vectors.

Theory basis:

- Earth-like worlds should not copy modern Earth, but their hypsometry,
  shelves, abyssal plains, trenches, ridge highs, land components, island size
  distribution, and coastline complexity should occupy plausible envelopes.

Implementation focus:

- Use equal-area statistics for global comparisons.
- Use distance-to-coast and distance-to-boundary profiles for shelves, slopes,
  abyssal plains, trenches, and nearshore sanity.
- Treat map PNGs as review evidence; the pass/fail logic belongs in metric
  fixtures.

Microbenchmark enrichment:

- Extend `P77.real_earth_hypsometry_extraction` with planform, shelf/slope,
  nearshore, island, and coastline metrics.
- Extend `P78.generated_hypsometry_envelope` with the same schema at 900,
  2500, and high-resolution audit scales.

Optimization targets:

- preserve current nearshore superdeep and far-ocean shoal safeguards;
- keep high terrain present but area-limited;
- prevent sea-level solving from exposing unsupported linear chains as
  mainland;
- preserve wet/dry variants as plausible alternatives rather than forcing all
  presets into modern-Earth land fraction.

### R2. Continental Province, Basement, Crust, and Sediment Evidence

Purpose:

- Resolve the main inland failure: broad continents should be composites of
  natural provinces, not single high smooth platforms.

Sources:

- USGS/NPS physiographic divisions and descriptions.
- Global geologic province and tectonic plate compilations.
- CRUST1.0 crustal thickness and sediment structure.
- GLiM lithology classes.
- NOAA total sediment thickness for margin and ocean checks.

Theory basis:

- Shields and cratons are old and stable but not automatically high plateaus.
- Covered platforms, intracratonic basins, foreland basins, failed rifts, and
  passive-margin lowlands are low because of accommodation, sediment, thermal
  history, and basement structure.
- Old sutures and mobile belts can remain natural boundaries after active
  mountain building ends.

Implementation focus:

- Promote province IDs/classes and parent process IDs into production fields
  rather than treating `terrain.continental_detail` as the only source of
  truth.
- Make base elevation consume province template, crust thickness, basement
  age, lithology/sediment proxies, and parent process before final noise.
- Keep high platforms constrained by explicit plateau/orogen/plume support.

Microbenchmark enrichment:

- Extend `P79.province_reference_graph_extraction` and
  `P80.generated_major_continent_province_graph`.
- Continue using `P83.crust_sediment_province_coupling` as the synthetic
  process-ordering reference.

Optimization targets:

- every large continent has several province classes;
- largest internal province share is capped for large continents;
- lowlands and basins exist on every large continent unless a documented
  special variant excludes them;
- high flat interiors require explicit plateau support and pass area caps.

### R3. Plate Boundary, Wilson-Cycle, and Deep-Time Evidence

Purpose:

- Ensure terrain-driving objects are geologically plausible across Earth
  history, not just present-day boundaries.

Sources:

- GPlates and EarthByte reconstruction workflows/data models.
- PB2002 present-day plate boundaries.
- World Stress Map.
- GEM Global Active Faults.
- Existing Aevum archive frames and Wilson-cycle objects.

Theory basis:

- A single Earth-like generator must support early hot-lithosphere variants,
  craton assembly, supercontinent cycles, ocean opening and closing, passive
  margin inheritance, subduction arcs, sutures, and old-orogen decay.
- Present-day absent features can still be required because they existed in
  the history that creates current continents.

Implementation focus:

- Make boundary process objects and lineage IDs the parents for terrain
  provinces.
- Preserve spherical continuity across the antimeridian in object logic, not
  only in raster rendering.
- Keep archive continuity for basin, margin, gateway, suture, orogen, and
  terrane identities.

Microbenchmark enrichment:

- Extend `P81.boundary_process_geometry_reference` for transform/ridge
  segmentation, fracture-zone adjacency, diffuse deformation belts, and trench
  polarity.
- Extend `P82.wilson_cycle_lifecycle_reference` for inherited passive margins,
  sutures, old orogens, and reactivated rifts.

Optimization targets:

- no major continent, seaway, or mountain belt should be chosen only by random
  seed texture;
- random seed may perturb geometry and texture, but parent events and process
  classes must be explainable by object state;
- archive fields should change only when a lifecycle transition justifies it.

### R4. Mountains, Plateaus, Rifts, Margins, and Special Landform Evidence

Purpose:

- Cover the real-Earth feature families that make continents read as natural:
  active margins, collision belts, plateaus, old subdued orogens, rift systems,
  passive-margin escarpments, volcanic/LIP provinces, hotspot tracks, arcs,
  back-arcs, marginal seas, and gateways.

Sources:

- GMBA mountain inventory.
- Global geologic province compilations.
- ETOPO/GEBCO relief and bathymetry profiles.
- Plate-boundary and reconstruction sources from R3.

Theory basis:

- Mountains can be narrow active ranges, broad collision plateaus, old eroded
  belts, volcanic arcs, rift-shoulder ranges, or extensional ranges.
- Plateaus require support and finite lifecycle: collision thickening, dynamic
  support, volcanic/LIP construction, or inherited highstanding crust.
- Rift-to-passive-margin systems should form ordered sequences of rift axis,
  shoulders, escarpment, coastal lowland, shelf, slope, rise, and abyssal
  plain.

Implementation focus:

- Promote landform objects from classification products into terrain drivers.
- Add parent-process requirements, area caps, age/decay fields, and lineage IDs
  for high-relief or special landforms.

Microbenchmark enrichment:

- Continue from `P87.mountain_inventory_expression`.
- Continue from `P88.rift_margin_escarpment_sequence`.
- Implement `P89.plateau_area_cap_and_decay`.
- Add successors only after P89 records which plateau fields are missing in the
  current generator.

Optimization targets:

- highlands and mountains remain visible but not line-noise dominated;
- plateaus are present, parented, finite, and decaying;
- passive-margin lowlands are large enough to matter but tied to shelf/sediment
  evidence;
- rift shoulders and escarpments are region objects rather than checkerboard
  cells.

### R5. Drainage, Erosion, and Source-to-Sink Evidence

Purpose:

- Move surface expression from isolated cell responses to drainage basin,
  erosion, and sediment-routing behavior.

Sources:

- HydroSHEDS, HydroBASINS, and HydroRIVERS.
- NOAA sediment thickness.
- GLiM lithology.
- ETOPO relief profiles.

Theory basis:

- Drainage divides should generally follow highlands, shields, old sutures,
  rift shoulders, and escarpments.
- Erosion lowers old orogens and mature interiors while preserving inherited
  province boundaries.
- Sediment should move from mountain/highland sources to foreland,
  intracratonic, passive-margin, delta/fan, shelf, and abyssal sinks.

Implementation focus:

- Build provisional drainage and sediment paths from province/base-relief
  surfaces before fine terrain texture.
- Keep climate-dependent discharge as a placeholder until climate and monsoon
  redesign resumes.

Microbenchmark enrichment:

- Continue from `P84.source_to_sink_sediment_budget`.
- Continue from `P85.drainage_divide_province_alignment`.
- Continue from `P86.old_orogen_erosion_decay`.

Optimization targets:

- sediment budget conservation within tolerance;
- drainage basins contiguous at region scale;
- old orogens lose relief while retaining boundary memory;
- sediment sinks lower basins and margins without turning all interiors into
  flat blankets.

### R6. Real-Earth Case-Study Evidence Packets

Purpose:

- Use named Earth regions as feature-class exemplars while avoiding exact map
  replay.

Initial case-study packets:

- North America: shield, covered platform, Cordilleran active margin,
  Appalachians as old orogen, foreland and intracratonic basins, Atlantic
  passive margin.
- South America: Andes active margin, foreland basin, craton/platform
  interior, Atlantic passive margin.
- Africa: craton/platform interiors, East African Rift, passive margins,
  volcanic/LIP highlands, sedimentary basins.
- Eurasia: Alpine-Himalayan collision belt, Tibetan-style plateau, old sutures,
  large platforms and basins.
- India and Australia: old shields/platforms, passive margins, rifted margins,
  Deccan-style volcanic plateau context.
- Ocean basins: Atlantic passive margins and ridge system, Pacific trenches
  and arcs, Indian Ocean ridges and fracture zones.

For each packet, archive:

- feature classes present;
- parent processes and lifecycle states;
- expected adjacency graph;
- elevation/relief/sediment profile sketches;
- current Aevum equivalent and residual gap;
- microbenchmark candidate or existing benchmark mapping.

Exit criteria:

- Every required Earth feature family has at least one case-study packet and at
  least one deterministic test fixture or generated-world audit.

### R7. Integrated Calibration and Promotion Matrix

Purpose:

- Combine evidence packets into a promotion decision without tuning one metric
  at the expense of another.

Implementation focus:

- Run small deterministic fixtures first, then generated-world audits, then
  8000-cell visual/contact-sheet review.
- Use parameter sweeps only after the theory owner and target metric are
  explicit.
- Record the parameter touched, expected response, observed response, and
  regression risks for every tuning pass.

Candidate benchmark:

- `P91.integrated_real_earth_morphology_promotion_audit`.

Promotion hierarchy:

1. Process/object correctness.
2. Archive and lifecycle continuity.
3. Province/crust/sediment/drainage ordering.
4. Hypsometry, bathymetry, and planform envelopes.
5. Compiler/render agreement.
6. PNG/contact-sheet visual reasonableness.

Exit criteria:

- Modern Earth-like worlds express the required real-Earth feature classes for
  geological reasons.
- Wet, dry, ocean, and continental variants remain intentionally different but
  do not violate core process logic.
- Remaining failures are named residuals with owners, not ambiguous visual
  objections.

## Theory Modules To Archive During Research

Each theory module should become a short design note before production code is
changed.  The note should identify parent processes, required fields, lifecycle
states, expected adjacency, terrain response, sediment response, and benchmark
acceptance.

1. **Stable continents and basement memory**
   - cratonization, shields, covered platforms, basement rework, lithosphere
     thickness, low strain, old age, and why stable does not mean uniformly
     high.
2. **Interior basins and platform topography**
   - sag basins, failed rifts, sediment cover, dynamic warping, accommodation,
     low relief, and broad lowland modes.
3. **Orogens and plateaus**
   - active collision, crustal thickening, foreland flexure, old-orogen decay,
     plateau support, erosion, and area caps.
4. **Rifts and passive margins**
   - extensional stress, inherited weakness, rift shoulders, aulacogens,
     thermal subsidence, passive-margin coastal plains, shelves, slopes, and
     sediment wedges.
5. **Active margins and ocean basins**
   - trench polarity, forearc/accretionary prism, volcanic arc, back-arc basin,
     ridge-transform-fracture fabric, abyssal plains, and age-derived
     bathymetry.
6. **Drainage and sediment**
   - drainage divides, basin hierarchy, erosion rates, transport capacity,
     source-to-sink conservation, delta/fan formation, and climate placeholders.
7. **Plumes, LIPs, and dynamic topography**
   - mantle heat anomaly, plume tracks, LIP finite area and decay, dynamic
     swells, and distinction from collision plateaus.
8. **Deep-time coverage**
   - Hadean/Archean hot lithosphere, proto-cratons, greenstone-like belts,
     Proterozoic supercontinents, Phanerozoic Wilson cycles, and present-day
     mature plate tectonics.

## Implementation Roadmap

### Stage A. Reference Data Preparation

Deliverables:

- source ledger schema and first source inventory;
- scripts or diagnostic functions for extracting small derived metrics;
- downsampled or summarized reference metric fixtures;
- no production behavior changes.

Candidate suites:

- `P76.reference_source_ledger_schema`
- `P77.real_earth_hypsometry_extraction`
- `P79.province_reference_graph_extraction`

Exit criteria:

- reference metrics can be reproduced offline from stored small fixtures;
- licenses and large-data storage policy are explicit;
- generated-world metrics use the same schema as reference metrics.

### Stage B. Generated-World Diagnostic Tightening

Deliverables:

- generated-world gates that compare against extracted reference envelopes;
- per-continent province graph metrics;
- richer failure categories: missing process, wrong adjacency, wrong amplitude,
  wrong scale, wrong lifecycle, renderer/compiler mismatch.

Candidate suites:

- `P78.generated_hypsometry_envelope`
- `P80.generated_major_continent_province_graph`
- `P81.boundary_process_geometry_reference`

Exit criteria:

- failures point to a process layer rather than a vague "map looks wrong";
- current Earth-like runs identify which province/terrain process is blocking
  promotion.

### Stage C. Province Architecture Production Layer

Deliverables:

- first-class production province graph per large continent;
- stable `tectonics.continental_province_id/code` or equivalent fields;
- parent process IDs for shield, platform, basin, old orogen, active orogen,
  foreland, rift, passive-margin lowland, and volcanic plateau objects;
- terrain base elevation derived from province templates before final detail.

Candidate suites:

- `P82.wilson_cycle_lifecycle_reference`
- `P83.crust_sediment_province_coupling`
- `P87.mountain_inventory_expression`
- promoted successors for P72/P74.

Exit criteria:

- large continents show multiple natural provinces;
- inland high/flat platform artifacts are reduced for process reasons;
- highlands and plateaus are area-limited and parented;
- passive margins, rifts, forelands, and old sutures remain recognizable.

### Stage D. Drainage, Erosion, and Sediment Coupling

Deliverables:

- drainage divide and basin generation from province/base relief;
- conservative sediment source-to-sink budget;
- foreland, intracratonic, passive-margin, delta/fan, and abyssal sinks;
- erosion decay for old orogens and surface lowering for stable platforms.

Candidate suites:

- `P84.source_to_sink_sediment_budget`
- `P85.drainage_divide_province_alignment`
- `P86.old_orogen_erosion_decay`
- updated E2 delta/fan and E5 glacial/surface-process guards.

Exit criteria:

- drainage and sediment affect terrain without land-mask regression;
- basins are not uniform sediment blankets;
- old orogens remain boundaries while relief decays.

### Stage E. High-Resolution Release and Promotion Audit

Deliverables:

- 8000/24000-cell Earth-like audit with the new reference metrics;
- optional 72000-cell deployment review only after profiling;
- contact sheets for elevation, province IDs, province classes, landform
  objects, drainage, sediment, crust age, ocean provinces, and compiled hex map.

Candidate suites:

- P75 is complete and should remain the promotion-decision baseline for the
  current P70-P74 evidence.
- Later promotion suites should replace weak legacy gates only when their
  province-architecture successors are stronger.

Exit criteria:

- modern Earth-like reference is inside broad feature-class envelopes;
- wet/dry variants remain plausible variants, not forced Earth copies;
- old local cleanup passes are no longer necessary for main correctness claims.

## Optimization Targets

Initial targets are deliberately envelopes, not exact constants.  Tighten them
after Stage A extraction.

Global/topographic targets:

- land fraction near the Earth-like envelope used by current release gates;
- multiple major land components in modern Earth-like mode;
- largest land component below previous P69/P21 blocker values;
- land ribbon fraction and narrow-neck density below current high-resolution
  blocker levels;
- land elevation mean and p95 within broad Earth-like ranges without uniform
  interior uplift;
- positive lowland share below 500 m and 1000 m on every large continent;
- highlands above 2500 m present but area-limited and process-parented;
- ocean shelf/slope/abyss/trench/ridge fractions preserved from current
  bathymetry sanity gates.

Province targets:

- at least three meaningful province classes per large continent;
- largest internal province fraction capped on large continents;
- shield/platform/basin/orogen/rift/passive-margin/foreland classes visible
  where their parent processes exist;
- province boundaries stable across archive frames unless a lifecycle event
  legitimately overprints them;
- no major province created only by random texture or final elevation
  thresholding.

Terrain and compiler targets:

- highlands, mountains, and plateaus must overlap parent process objects;
- basins and passive-margin lowlands must be low for accommodation reasons;
- drainage divides should align with province boundaries or highland objects;
- hex compiler land/ocean, shelf/deep-ocean, and elevation sign consistency
  must stay at current hard-gate quality;
- rendered PNGs are audit evidence, not sole acceptance.

## P101+ Real-Earth Comparison Repair Archive

Archived: 2026-06-27

This is the detailed plan for the next work after P100.  The priority is still
plate/land/terrain realism, not climate, ocean circulation, or monsoon coupling.
The plan turns real Earth landforms into source packets, theory notes,
microbenchmarks, implementation owners, and optimization targets before any new
large production rewrite.

### Current Residual Baseline

Use the after-P99/P100 evidence as the frozen starting point:

- `P100.integrated_reaudit_and_promotion_gate` passes and verifies all P93-P99
  repair suites.
- Default Earth-like promotion remains blocked by four named blockers:
  `p69_earthlike_reference_needs_calibration`,
  `p90_current_world_residuals_unresolved`,
  `crust_sediment_residuals_unresolved`, and
  `planform_residuals_unresolved`.
- P90 after-P99 records `19` gaps: `10` planform, `1` crust/sediment, and `8`
  asset-review entries.
- The non-asset residuals are concentrated in:
  land fraction too low (`0.2311`), largest land component too dominant
  (`0.7837`), major land component count too low (`2`), land ribbon fraction too
  high (`0.5817`), largest-component coastline complexity slightly too high
  (`8.6159`), land elevation p95 slightly too low (`1471.8 m`), and high-flat
  interior share too high (`0.0625` of continental land).

P101 must not clear these by relaxing reference envelopes.  It should either
repair production behavior or record a narrower blocker with a stronger theory
justification.

### Archive Packet Template

Every phase below should create or update an evidence packet.  A packet is small
and commit-friendly; raw source rasters/vectors remain outside git.

Required fields:

- `packet_id` and owner stage.
- source ids, version notes, license notes, and local raw-data policy.
- extraction or curation method, including projection and equal-area weighting.
- theory claims defended by the packet.
- derived metrics with units and expected ranges.
- synthetic fixture and generated-world audit using the same schema.
- optimization targets and forbidden tuning shortcuts.
- expected PNG/contact-sheet evidence.
- residual policy: promotion blocker, warning, or deferred future work.

### Phase 0. Baseline Reproduction and Failure Attribution

Goal:

- Reproduce P90/P91/P100 residuals from a clean current run.
- Attach each residual to a code owner: continent assembly, land exposure,
  sea-level solve, province/crust/sediment coupling, terrain base elevation, or
  compiler/render.

Theory basis:

- A visual landform defect is actionable only after it is expressed as a process
  failure or a measurement failure.
- Planform and inland-elevation failures should be separated before tuning.

Implementation method:

- Add `P101.planform_crust_sediment_residual_repair` with a read-only baseline
  sub-benchmark first.
- Store current residual metrics in the P101 summary, including generated-world
  globals that reveal which repair systems fired.
- Keep this phase production-read-only unless the baseline cannot be reproduced.

Test details:

- Generated worlds: Earth-like at `900` cells for CI and `2500` cells for
  intermediate confidence.
- Required metrics: P90 gap ids, owner counts, category counts, direct planform
  metrics, high-flat interior share, per-component area/lowland/high-flat
  summary, and compiler consistency.
- Acceptance: metrics match the archived after-P99/P100 residual shape and every
  residual has one primary owner plus a regression-risk note.

Optimization target:

- No optimization yet; this phase freezes the failing surface.

### Phase 1. Reference Evidence Expansion

Goal:

- Enrich real-Earth reference packets so P101 and later phases tune toward
  geological envelopes, not subjective map appearance.

Source groups:

- ETOPO/GEBCO/Natural Earth for global hypsometry, bathymetry, coastlines,
  islands, shelves, slopes, abyssal plains, ridges, trenches, and component
  geometry.
- USGS/NPS/global tectonic-province sources for continent-internal province
  graphs.
- CRUST1.0, GLiM, and global sediment-thickness sources for basement, lithology,
  crust thickness, and sediment accommodation.
- GPlates/EarthByte/PB2002/World Stress Map/GEM active faults for boundary
  process classes and deep-time case studies.
- HydroSHEDS/HydroBASINS/HydroRIVERS and GMBA for drainage, erosion, mountain
  inventory, and source-to-sink context.

Theory basis:

- Modern Earth is a calibration set for feature classes and distributions, not a
  literal target map.
- First-order geography comes from coupled planform, crustal support, sediment,
  boundary lifecycle, and surface-process history.

Implementation method:

- Extend the P76 ledger and P77/P79 reference extractors before changing
  production coefficients.
- Store small JSON fixtures with checksums and source ids.

Test details:

- Reference track: extraction reproducibility, source id stability, checksum
  stability, and expected broad values.
- Generated-world track: same schema at 900/2500/high-resolution audit scales.

Optimization targets:

- Land fraction and component envelopes remain broad Earth-like ranges.
- Island/ribbon/coastline metrics must be resolution-aware.
- Shelf/slope/abyss/trench/ridge metrics must distinguish nearshore depth errors
  from valid subduction trenches and mid-ocean ridge highs.

### Phase 2. Planform Mechanism Repair

Goal:

- Replace the remaining long-ribbon / two-landmass / overdominant-continent
  behavior with process-backed continent architecture.

Theory basis:

- Linear geological features are valid as arcs, rifts, sutures, margins, and
  volcanic chains; they become invalid when they are the exposed continental mask
  without broad continental interiors or component-scale support.
- Earth-like worlds need several major exposed components or a documented
  supercontinent variant.

Implementation method:

- Audit continent assembly and late land exposure together: initial nuclei,
  later merging, continental conservation, sea-level solve, shelf/platform
  exposure, island-arc exposure, and final coastline smoothing.
- Prefer object-supported broad interiors, embayments, shelves, and passive
  margins over arbitrary ocean-to-land flips.
- Preserve spherical continuity and avoid antimeridian-special logic.

Candidate P101 sub-benchmarks:

- `P101.planform_residual_baseline`
  - metrics: current P90 planform gaps and per-component summaries;
  - acceptance: baseline reproduced and owner-attributed.
- `P101.component_architecture_repair`
  - metrics: land fraction, major component count, largest component share,
    component lowland/high-flat shares;
  - acceptance: land fraction enters the Earth-like envelope without merging into
    a single dominant continent.
- `P101.ribbon_coastline_repair`
  - metrics: land ribbon fraction, continental-crust ribbon fraction, largest
    coastline complexity, island area distribution;
  - acceptance: exposed mainland ribbons fall below the current blocker while
    valid island arcs and rift margins remain parented.

Optimization targets:

- land fraction in the broad Earth-like envelope, currently `0.25-0.33`;
- major component count in the broad Earth-like envelope, currently `4-14`;
- largest land component share below the current blocker envelope upper bound,
  currently `0.60`;
- land ribbon fraction below the current blocker envelope upper bound, currently
  `0.35`;
- largest-component coastline complexity no greater than the current envelope
  upper bound, currently `8.0`, unless the reference extraction revises it.

### Phase 3. Crust, Sediment, and Interior Elevation Repair

Goal:

- Fix high, flat continental interiors without deleting valid highlands,
  mountains, and plateaus.

Theory basis:

- Stable cratons and shields preserve old basement and low strain; they are not
  automatically high tablelands.
- Covered platforms, intracratonic basins, foreland basins, failed rifts, and
  passive-margin lowlands lower or partition interiors through accommodation,
  sediment, thermal history, and inherited basement structure.
- High plateaus require explicit support: collision thickening, dynamic support,
  LIP/plume construction, or inherited highstanding crust with lifecycle limits.

Implementation method:

- Make final base elevation consume province class, crust thickness, sediment
  accommodation, basement age/rework, drainage/erosion, and parent process before
  final texture.
- Preserve lowland and basin expression per large continent.
- Add relief or lowering only through process owners; do not add noise solely to
  break a metric.

Candidate P101/P102 sub-benchmarks:

- `P101.crust_sediment_high_flat_repair`
  - metrics: high-flat interior fraction by world, continent, and province;
  - acceptance: high-flat blocker clears while supported plateaus and orogens
    remain visible.
- `P102.interior_lowland_basin_ordering`
  - metrics: shield/platform/basin/old-orogen/passive-margin lowland elevation,
    local relief, crust thickness, sediment thickness, and parent object
    coverage;
  - acceptance: basins and passive margins are lower than shields/platforms,
    active orogens are higher and rougher, and covered platforms are not uniform
    elevated slabs.

Optimization targets:

- high-flat interior share below the current P90 blocker threshold, currently
  `0.02` of continental land;
- positive lowland share below `500 m` and `1000 m` on each large continent;
- land p95 high enough to preserve mountains, currently at least `1500 m` in the
  broad envelope;
- no highland or plateau area without active, old, plume/LIP, rift-shoulder, or
  inherited-crust support.

### Phase 4. Province Graph and Natural Continental Boundaries

Goal:

- Ensure each major continent reads as a composite of natural provinces with
  meaningful internal boundaries.

Theory basis:

- Real continents contain shields, platforms, basins, old orogens, active
  margins, rift systems, forelands, passive margins, volcanic provinces, and
  inherited sutures in characteristic adjacency graphs.
- Natural boundaries can remain visible after active tectonics fades.

Implementation method:

- Treat `tectonics.continental_provinces` and production province objects as the
  terrain driver, not a late classification product.
- Preserve province boundaries through archive lineage ids and parent process
  references.
- Link drainage divides and sediment sinks to province boundaries only after the
  base province graph is stable.

Test details:

- Extend P80/P94/P96 generated-world audits with per-major-continent province
  count, largest province share, adjacency edge coverage, parent process
  coverage, and boundary persistence.
- Add visual review assets for province id/class maps and elevation overlays.

Optimization targets:

- every large continent has several province classes;
- no single platform/shield province dominates a large continent unless the
  preset explicitly asks for a small or ancient craton-dominated world;
- internal boundaries are visible through relief, basin, drainage, lithology, or
  sediment differences without becoming random checkerboard noise.

### Phase 5. Landform Family and Case-Study Coverage

Goal:

- Ensure real Earth feature families can be generated through model mechanisms.

Case-study packets:

- North America-style package: shield, covered platform, Cordilleran active
  margin, Appalachian old orogen, foreland/intracratonic basins, Atlantic passive
  margin.
- South America-style package: Andes active margin, foreland basin,
  craton/platform interior, Atlantic passive margin.
- Africa-style package: cratons/platforms, East African Rift, passive margins,
  volcanic highlands, sedimentary basins.
- Eurasia-style package: collision belt, Tibetan-style plateau, old sutures,
  large platforms and basins.
- India/Australia-style package: old shields, passive margins, rifted margins,
  Deccan-style volcanic context.
- Ocean-basin package: Atlantic ridge/passive margins, Pacific trenches/arcs,
  Indian Ocean ridges/fracture zones.

Test details:

- Each case-study packet maps to at least one synthetic fixture and one
  generated-world audit.
- Required measurements: feature classes, parent processes, adjacency graph,
  elevation/relief/sediment profile, lifecycle state, and current Aevum residual.

Optimization targets:

- all required feature families in `docs/EARTH_GEOMORPHOLOGY_COVERAGE.md` are
  covered either by generated-world expression or a named fixture with a
  production path;
- fixture-only features remain blockers until they appear in ordinary
  generated-world audits or are intentionally gated to a preset variant.

### Phase 6. Multi-Resolution Calibration and Asset Review

Goal:

- Tune with numeric evidence and inspect PNGs only after the metrics indicate the
  same failure mode has moved.

Implementation method:

- Run every major change at 900 cells first, then 2500 cells, then 8000-cell
  asset review.
- Use high-resolution assets to detect visual artifacts that metrics may miss:
  long straight chains, uniform high slabs, speckled archipelagos, nearshore
  depth inversions, compiler/elevation mismatches, and province-map/elevation-map
  disagreement.

Required asset set:

- `elevation.png`
- `terrain_provinces.png`
- `continental_detail_provinces.png`
- `ocean_depth_provinces.png`
- `plates.png`
- `crust_age.png`
- `history.png`
- `timeline.png`
- `hexmap.png`
- P91/P101 contact sheets

Test details:

- Asset existence is not enough.  P101+ summaries should record which assets are
  for review, which residuals they correspond to, and whether compiler
  consistency still passes.

Optimization targets:

- no regression in cleared P93-P99 owner blockers;
- no compiler/render land-ocean sign mismatch;
- metric movement must match the intended process owner.

### Phase 7. Promotion Decision and Next-System Boundary

Goal:

- Decide whether plate/land/terrain can move from repair mode to downstream
  systems.

Implementation method:

- After P101/P102-style repairs, rerun P90, P91, and P100.
- Promote only if planform and crust/sediment residual blockers clear without
  reintroducing province graph, boundary lifecycle, drainage/erosion, landform,
  or bathymetry/margin blockers.
- Keep climate, ocean-current, and monsoon redesign paused until this gate is
  green or the remaining failures are explicitly outside plate/terrain ownership.

Exit criteria:

- root non-asset P90 gaps for planform and crust/sediment are cleared or
  replaced by narrower accepted residuals;
- P91 reports no unresolved root owner blockers other than the separate P69
  calibration task;
- P100 or successor promotion gate records `release_gate_allowed: True`, or a
  precise remaining blocker set with the next owner packet.

## P103+ Source-Corpus Enrichment and Repair Planning Archive

Archived: 2026-06-27

This archive expands the P101+ repair sequence with the requested research
plan: collect substantial real-Earth evidence in stages, then use it to enrich
test details, theory basis, implementation methods, and optimization targets.
It does not replace P101/P102.  It makes the later P103+ repair packets harder
to game by tying every metric and production change to a named evidence source
and a geological mechanism.

The immediate engineering priority remains real-Earth landform comparison and
plate/terrain repair.  Climate, ocean circulation, and monsoon redesign stay
paused except for regression fixes.

### External Source Corpus

Raw source rasters and vectors must remain outside git.  The repository should
store source ids, version/access metadata, checksums for derived fixtures, small
JSON metric fixtures, and extraction scripts only.

Tier 1 sources for global shape and relief:

- `NOAA_ETOPO_2022`: NOAA NCEI ETOPO 2022 Global Relief Model.
  URL: `https://www.ncei.noaa.gov/products/etopo-global-relief-model`
  Use: land/ocean hypsometry, topographic tails, shelf/slope/abyss profiles,
  nearshore depth sanity, and global elevation bins.
- `GEBCO_2026_GRID`: GEBCO gridded bathymetry.
  URL: `https://www.gebco.net/data-products/gridded-bathymetry-data`
  Use: ocean floor depth provinces, trenches, ridges, abyssal plains, and
  under-ice/land-ocean terrain cross-checks.
- `NATURAL_EARTH_10M_PHYSICAL`: Natural Earth physical vectors.
  URL: `https://www.naturalearthdata.com/downloads/10m-physical-vectors/`
  Use: coastline, land/ocean polygons, island/component geometry, coastline
  complexity, and planform generalization checks.

Tier 2 sources for crust, sediment, lithology, and provinces:

- `CRUST1_0`: CRUST1.0 global crustal model.
  URL: `https://ds.iris.edu/ds/products/emc-crust10/`
  Use: crust thickness, sediment thickness baseline, continent/ocean crust
  ordering, and high-standing crust plausibility.
- `ECM1`: Earth Crustal Model 1.
  URL: `https://www.earthcrustmodel1.com/`
  Use: cross-check crustal type, sediment, crystalline crust thickness, and
  continental support ranges.
- `GLIM_GLOBAL_LITHOLOGY`: Global Lithological Map.
  URL: `https://www.geo.uni-hamburg.de/en/geologie/forschung/aquatische-geochemie/glim.html`
  Use: lithology, shield/platform/basin surface expression, and erosion
  resistance context.
- `NOAA_TOTAL_SEDIMENT_THICKNESS`: NOAA total sediment thickness for oceans and
  marginal seas.
  URL: `https://www.ncei.noaa.gov/products/total-sediment-thickness-oceans-seas`
  Use: passive-margin sediment wedges, basin fill, shelf/rise profiles, and
  source-to-sink acceptance ranges.
- `USGS_WORLD_GEOLOGIC_MAPS`: USGS world geologic maps.
  URL: `https://www.usgs.gov/centers/central-energy-resources-science-center/science/world-geologic-maps`
  Use: broad lithologic/province context and case-study checks where global
  province graph data are incomplete.

Tier 3 sources for plate boundaries, active deformation, and deep-time context:

- `GPLATES_EARTHBYTE_MODELS`: EarthByte/GPlates global and regional plate
  motion models.
  URL: `https://www.earthbyte.org/category/resources/data-models/global-regional-plate-motion-models/`
  Use: deep-time plate reconstruction, Wilson-cycle examples, ridge/trench
  lifecycles, sutures, and old-orogen inheritance.
- `GPLATES_WEB_SERVICE_MODELS`: GPlates Web Service model catalog.
  URL: `https://gwsdoc.gplates.org/models/`
  Use: reproducible model metadata and optional automated reconstruction checks.
- `PB2002_PLATE_BOUNDARIES`: Bird PB2002 plate boundary model.
  URL: `https://peterbird.name/publications/2003_pb2002/2003_pb2002.htm`
  Use: modern ridge/transform/trench/collision boundary classes and boundary
  length envelopes.
- `GEM_GLOBAL_ACTIVE_FAULTS`: GEM Global Active Faults database.
  URL: `https://www.globalquakemodel.org/product/active-faults-database`
  Use: active fault trace geometry, kinematics, slip-rate metadata, and diffuse
  deformation case-study context.
- `WORLD_STRESS_MAP`: World Stress Map.
  URL: `https://www.world-stress-map.org/`
  Use: present-day stress regime context, intraplate deformation plausibility,
  and active-margin/diffuse-zone validation.

Tier 4 sources for drainage, erosion, mountains, and source-to-sink structure:

- `HYDROBASINS`: HydroBASINS / HydroSHEDS catchment hierarchy.
  URL: `https://www.hydrosheds.org/products/hydrobasins`
  Use: drainage basin hierarchy, divide-to-outlet geometry, basin scale
  distributions, and continent-internal drainage partitioning.
- `HYDRORIVERS`: HydroRIVERS.
  URL: `https://www.hydrosheds.org/products/hydrorivers`
  Use: river reach density, outlet routing, and source-to-sink graph examples.
- `MERIT_HYDRO`: MERIT Hydro global hydrography.
  URL: `https://global-hydrodynamics.github.io/MERIT_Hydro/`
  Use: flow direction, upstream area, HAND, river width, and drainage relief
  cross-checks where HydroSHEDS limitations matter.
- `GMBA_MOUNTAIN_INVENTORY`: GMBA Mountain Inventory.
  URL: `https://www.earthenv.org/mountains`
  Use: mountain range hierarchy, mountain area scale, range clustering, and
  case-study family coverage.

### Source Collection Phases

Phase S0. Ledger hardening:

- Extend P76/P102 source metadata with `source_url`, `access_date`,
  `version_or_release`, `license_status`, `raw_storage_policy`,
  `derived_fixture_path`, `checksum_status`, and `intended_metric_groups`.
- Add a source-ledger schema check that fails on missing version/access/license
  fields for any source used by a benchmark.
- Record sources that are literature-only or metadata-only rather than
  pretending raw GIS extraction has been done.

Phase S1. Global relief and planform extraction:

- Build or regenerate small fixtures for global hypsometry, ocean bathymetry,
  land component geometry, island area distribution, coastline complexity,
  shelf/slope/rise/abyss fractions, nearshore-to-offshore depth gradients, and
  antimeridian continuity.
- Use equal-area weighting and store projection/resolution assumptions in the
  fixture metadata.
- Produce a fixture-vs-generated comparison using the same metric keys as P90,
  P91, P101, and P102.

Phase S2. Crust, sediment, lithology, and province extraction:

- Extract coarse reference envelopes for crust thickness by province setting,
  sediment thickness by passive margin/foreland/intracratonic basin, exposed
  lithology by province class, and old basement/platform/basin relationships.
- Preserve separate rules for stable shield, covered platform, basin, old
  orogen, active orogen, plateau, rift, passive margin, and shelf/rise.
- Add a fixture proving the core rule: old stable crust can be low or moderate
  relief and must not automatically become a high flat tableland.

Phase S3. Boundary and Wilson-cycle extraction:

- Derive modern boundary length fractions and adjacency rules from PB2002/GEM
  sources, then cross-check deep-time lifecycle examples with GPlates/EarthByte
  models.
- Track ridge/transform/trench/collision/suture/passive-margin/rift/plume/LIP
  as objects with lifecycle state rather than final raster labels.
- Add antimeridian and spherical-continuity tests to every boundary fixture.

Phase S4. Surface-process and mountain-family extraction:

- Use HydroBASINS/HydroRIVERS/MERIT Hydro to define drainage basin hierarchy,
  divide persistence, outlet routing, endorheic/through-flow distinctions, and
  basin-size envelopes.
- Use GMBA mountain polygons to define mountain range area, clustering, and
  hierarchy envelopes.  These are range-family targets, not exact modern Earth
  coordinate targets.
- Connect mountain and drainage fixtures to parent processes: collision,
  subduction, old-orogen relict, rift shoulder, plume/LIP, and inherited
  basement edge.

Phase S5. Case-study packet expansion:

- Maintain six continent/ocean packages: North America, South America, Africa,
  Eurasia, India/Australia, and ocean basins.
- Each package must list: feature classes, parent processes, adjacency graph,
  expected elevation/relief/sediment signals, lifecycle state, source ids,
  generated-world audit target, and current residual risk.
- Case studies are for feature-family coverage.  They must not force Aevum to
  clone modern Earth geography.

Phase S6. Production repair mapping:

- Before each production change, write the target benchmark name, owner layer,
  affected source packets, affected code entry points, forbidden shortcuts, and
  expected metric movement.
- After each production change, rerun the owner packet plus P101/P102 and the
  relevant cleared-owner regression suites.
- If a metric improves by changing the wrong process owner, keep the blocker
  open and record the mismatch.

Phase S7. Promotion and residual audit:

- Rerun the current-world gap inventory, integrated promotion audit, and
  rendered asset review after every major source-informed repair phase.
- Promotion requires that planform and crust/sediment blockers clear without
  reintroducing province graph, boundary lifecycle, drainage/erosion, landform
  expression, bathymetry/margin, or compiler/render blockers.
- If residuals remain, record their owner, source-packet status, and next
  repair packet.  Do not hide residuals behind wider envelopes unless the
  reference extraction itself justifies a threshold change.

### Enriched Microbenchmark Plan

The following P103+ benchmark sequence is reserved for execution after P102.
Names can be split if implementation risk requires smaller patches.

- `P103.planform_mechanism_repair`
  - Evidence: S1 plus P101/P102.
  - Code owners: continent assembly, land exposure, sea-level solve,
    component preservation, ribbon pruning, coastline smoothing.
  - Metrics: land fraction, major component count, largest component share,
    ribbon fraction, island area distribution, coastline complexity,
    antimeridian continuity, and compiler land/ocean agreement.
  - Target: move the generated Earth-like planform into the broad reference
    envelope without creating a single dominant continent or speckled islands.
- `P104A.continental_mosaic_object_expression`
  - Evidence: S2/S5 plus P79/P80/P94/P102 and generated-world map review.
  - Code owners: production continental province graph,
    `terrain.continental_detail`, terrain/elevation response, map compiler
    sampling, and province object diagnostics.
  - Metrics: per-major-continent province class count, largest internal
    province share, shield/platform/basin/rift/orogen/plateau coverage,
    province-code-to-detail agreement, relief variance inside major continents,
    object parent-process coverage, and rendered
    `continental_detail_provinces.png`/`elevation.png` readability.
  - Target: each large continent is a mosaic of persistent cratons, platforms,
    basins, old orogens, failed rifts, sutures, passive margins, and plateaus.
    Major province existence must come from parent process objects or inherited
    crust state, not from random texture or last-step recoloring.
- `P104B.crust_sediment_interior_elevation_repair`
  - Evidence: S2 plus P83/P96/P101.
  - Code owners: crust/sediment/province ordering, base elevation,
    isostatic support, erosion/lowland preservation.
  - Metrics: high-flat interior share, basin/platform/shield elevation
    ordering, sediment accommodation, crust-thickness support, land p95,
    lowland share per major continent, and unsupported highland area.
  - Target: remove high flat interiors while preserving the P104A continental
    mosaic and valid mountains, plateaus, shields, basins, and passive-margin
    lowlands.
- `P105.natural_province_boundary_expression`
  - Evidence: S2/S5 plus P79/P80/P94.
  - Code owners: production province graph, province lineage, terrain
    expression, drainage/sediment boundary coupling.
  - Metrics: province count per major continent, largest province share,
    class-edge coverage, boundary persistence, parent-process coverage, and
    non-checkerboard boundary visibility.
  - Target: every large continent contains multiple process-backed provinces
    separated by natural physiographic boundaries.
- `P106.boundary_and_wilson_cycle_consistency`
  - Evidence: S3 plus P81/P82/P95.
  - Code owners: boundary objects, Wilson-cycle archive, spreading centers,
    sutures, passive margins, rifts, trenches, and transforms.
  - Metrics: boundary process coverage, lifecycle transitions, object
    persistence, spherical continuity, and final terrain expression.
  - Target: generated modern geography is the surface expression of deep-time
    plate events rather than a late raster redraw.
- `P107.drainage_mountain_source_to_sink_coupling`
  - Evidence: S4 plus P84/P85/P86/P87/P97/P98.
  - Code owners: drainage basins, divides, erosion, sediment routing, mountain
    inventory, old-orogen decay, and landform lifecycle.
  - Metrics: basin hierarchy, divide/province alignment, mountain range area
    hierarchy, sediment budget closure, outlet/sink mapping, and eroded
    old-orogen relief.
  - Target: interior relief and lowlands are organized by drainage and
    sediment systems, not by per-cell noise.
- `P108.real_earth_case_study_family_gate`
  - Evidence: S5 plus all relevant prior packets.
  - Code owners: integration across planform, crust/sediment, boundary,
    drainage, landform, and bathymetry/margin layers.
  - Metrics: feature-family coverage by case-study package, generated-world
    occurrence, fixture-only gaps, parentless features, and asset review notes.
  - Target: the model can generate Earth-realistic feature families for the
    right process reasons, even when the map is not modern Earth.
- `P109.multi_resolution_asset_promotion_audit`
  - Evidence: S1-S5 plus P90/P91/P100 successor gates.
  - Code owners: diagnostics, renderer, compiler, release gate.
  - Metrics: 900/2500/8000-cell stability, metric drift by resolution, PNG
    completeness, visual residual linkage, compiler/elevation agreement, and
    promotion blocker count.
  - Target: decide whether plate/terrain can leave repair mode or whether the
    remaining blockers require another owner-specific packet.

### Theory Notes Required Before Production Tuning

Every new benchmark must include a short theory note with these fields:

- process claim: what geological mechanism should create the feature;
- expected geometry: broad component, belt, margin, basin, plateau, trench,
  shelf, fan, divide, or range hierarchy;
- lifecycle state: active, inherited, decaying, buried, drowned, uplifted, or
  reactivated;
- required parent object: plate boundary, province, crustal block, basin,
  drainage network, sediment sink, plume/LIP, or archive lineage;
- allowed randomness: local texture, perturbation, or ensemble variability only;
- forbidden shortcut: random major continent split, arbitrary sea-level dodge,
  visual-only recoloring, or threshold relaxation without source evidence;
- expected generated-world metric movement;
- regression risks to already-cleared owner layers.

### Optimization Discipline

- Optimize process parameters, not final rasters.  A repair that only flips
  cells to satisfy a component metric is invalid unless the flipped cells have a
  parent province, margin, shelf, basin, or boundary object.
- Random seeds may vary local texture and ensemble sampling.  They must not
  choose major landmasses, Wilson-cycle events, craton survival, province
  existence, or final promotion eligibility.
- Thresholds can change only after the reference fixture changes and records
  the source reason.
- Every owner repair must prove that it does not regress the already-cleared
  owner layers: province graph, boundary lifecycle, drainage/erosion,
  landform expression, bathymetry/margin, and compiler/render consistency.
- Rendered PNG review is required but subordinate to metrics.  A visual issue
  becomes actionable only after it maps to a metric, source packet, and owner
  layer.

## Progress Log

2026-06-26 - Plan archived

- Added this document as the detailed real-Earth geomorphology research,
  theory, implementation, optimization, and microbenchmark archive.
- Kept climate/ocean-current/monsoon redesign explicitly out of scope for this
  cycle.
- Anchored the next work to current P70-P74 assets and the pending P75
  promotion audit.
- Added candidate suites P76-P89 to convert source collection and theory work
  into reproducible diagnostics before production rewrites.

Next:

- Complete P75 release/promotion audit using the current P70-P74 evidence.
- Then begin Stage A with `P76.reference_source_ledger_schema` and
  `P77.real_earth_hypsometry_extraction`.

2026-06-27 - P75 release/promotion audit completed

- Added the executable `P75.release_and_promotion_audit` benchmark.
- P75 selects the archived passing P69 high-resolution audit, verifies complete
  8000-cell key assets, reruns P70-P74, and reruns P29/P48/P49/P68 legacy
  gates.
- P75 passes as a release audit but intentionally does not promote the current
  behavior to default Earth-like generation.
- Current blockers are:
  `reference_data_download_or_raster_extraction_pending`,
  `first_class_production_province_graph_pending`, and
  `p69_earthlike_reference_needs_calibration`.
- Archived output:
  `/Users/rayw/Projects/aevum/out_bench_p75_release_promotion_audit_20260627/`.

Next:

- Begin Stage A with `P76.reference_source_ledger_schema` and
  `P77.real_earth_hypsometry_extraction`.

2026-06-27 - P76 reference source ledger schema completed

- Added an offline source-ledger schema in
  `aevum/diagnostics/reference_source_ledger.py`.
- The ledger covers all P70 reference sources and records version notes,
  explicit license status, acquisition/extraction status, raw-data policy,
  local storage policy, projection/resolution notes, checksum status, and
  derived metric targets.
- The benchmark is deliberately metadata-only: external raw reference data is
  not downloaded or required by default.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p76_reference_source_ledger_schema_20260627/`.
- Current metrics:
  source count `18`, required field count `14`, missing required fields `0`,
  duplicate source IDs `0`, invalid phase references `0`.

Next:

- Implement `P77.real_earth_hypsometry_extraction`.

2026-06-27 - P77 real-Earth hypsometry extraction fixture completed

- Added a small derived hypsometry metric fixture at
  `data/reference/earth_hypsometry_fixture_20260627.json`.
- Added `aevum/diagnostics/real_earth_hypsometry.py` and the executable
  `P77.real_earth_hypsometry_extraction` benchmark.
- The fixture stores 12 area-weighted elevation/depth bins and derived metrics
  for land fraction, land-elevation mean/p95, highland tail, lowland shares,
  ocean mean depth, shelf/slope-rise/abyss fractions, and shelf-to-abyss depth
  separation.
- Raw ETOPO/GEBCO rasters are not committed; the fixture explicitly marks
  direct raster extraction and future regeneration as pending.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p77_real_earth_hypsometry_extraction_20260627/`.
- Current metrics:
  land fraction `0.292`, land mean elevation `835 m`, land p95 `2600 m`,
  high-land fraction above 2500 m `0.055`, shelf/slope-rise/abyss ocean
  fractions `0.056/0.117/0.636`, and envelope checks `8/8`.

Next:

- Implement `P78.generated_hypsometry_envelope`.

2026-06-27 - P78 generated hypsometry envelope completed

- Added generated-world hypsometry metrics and comparison against the P77
  fixture schema.
- The P78 gate uses a current 900-cell Earth-like world for current-code signal
  and archived P69 8000-cell evidence for high-resolution continuity.
- The generated world's core hypsometry now passes the broad reference
  envelope, but promotion is still blocked by recorded planform/trench residual
  metrics.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p78_generated_hypsometry_envelope_20260627/`.
- Current generated metrics:
  land fraction `0.240`, land mean elevation `859 m`, land p95 `1928 m`,
  high-land fraction above 2500 m `0.009`, shelf/slope-rise/abyss ocean
  fractions `0.091/0.190/0.308`, and shelf-to-abyss depth delta `3387 m`.

Next:

- Implement `P79.province_reference_graph_extraction`.

2026-06-27 - P79 province reference graph extraction completed

- Added a small derived real-Earth province graph fixture in
  `aevum/diagnostics/province_reference_graph.py`.
- The fixture is extracted from the P73 case-study calibration sketches and
  makes province class/process/adjacency expectations executable without
  committing raw GIS vectors.
- Added `P79.province_reference_graph_extraction` to the tectonics benchmark
  suite and pytest coverage for graph schema, coverage, connectivity, and raw
  vector policy.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p79_province_reference_graph_extraction_20260627/`.
- Current metrics:
  case count `5`, province nodes `29`, adjacency edges `25`, province classes
  `9`, parent processes `15`, source IDs `9`, class-edge types `16`, and
  missing required feature classes/processes/class edges `0/0/0`.
- Raw vector extraction remains a future regeneration step; P79 records it as
  pending rather than pretending this fixture is exact geography.

Next:

- Implement `P80.generated_province_graph_reference_comparison`.

2026-06-27 - P80 generated province graph reference comparison completed

- Added `aevum/diagnostics/generated_province_reference.py`.
- Added `P80.generated_province_graph_reference_comparison` to the tectonics
  benchmark suite.
- P80 maps generated detail classes and landform objects onto the P79 reference
  graph classes, then checks generated class coverage, mapped parent-process
  coverage, major-continent multi-class structure, dominant-class caps, and
  reference class-edge gaps.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p80_generated_province_graph_reference_comparison_20260627/`.
- Current metrics:
  generated/reference class count `8/9`, mapped/required parent processes
  `8/8`, generated/required class edges `23/9`, major continents `2`, minimum
  major-continent reference class count `7`, maximum largest reference-class
  fraction `0.473`, missing required feature classes/processes `0/0`, and
  unexpected missing reference classes/edges `0/0`.
- Recorded residuals:
  `volcanic_lip_plateau` and `rift_system|volcanic_lip_plateau`.
- Interpretation:
  P80 establishes a real comparison gate rather than another visual check.
  The generated world satisfies current required province/process coverage, but
  first-class production province IDs and LIP/volcanic plateau expression remain
  future work.

Next:

- Implement `P81.boundary_process_geometry_reference`.

2026-06-27 - P81 boundary process geometry reference completed

- Added `aevum/diagnostics/boundary_process_geometry.py`.
- Added `P81.boundary_process_geometry_reference` to the tectonics benchmark
  suite.
- P81 introduces an executable boundary-process geometry reference: source IDs,
  broad length-fraction envelopes, a deterministic synthetic spherical boundary
  network, and a current generated-world boundary comparison using the same
  schema.
- The synthetic fixture covers ridge, transform, subduction trench,
  collision/suture, diffuse deformation, passive margin, and continental rift
  process types.  It checks transform-offset adjacency, trench-active-margin
  adjacency, collision-diffuse adjacency, curved geometry, and antimeridian
  continuity.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p81_boundary_process_geometry_reference_20260627/`.
- Current metrics:
  source IDs `5`, synthetic process types `7`, synthetic length-envelope
  checks `7/7`, transform offset count `2`, antimeridian ridge components `1`,
  transform-near-ridge fraction `0.846`, trench-near-active-margin fraction
  `1.000`, generated-world process types `6`, generated boundary cells `103`,
  and generated boundary cell fraction `0.114`.
- Recorded residual:
  current generated Earth-like worlds still lack transform boundary masks.

Next:

- Implement `P82.wilson_cycle_lifecycle_reference`.

2026-06-27 - P82 Wilson-cycle lifecycle reference completed

- Added `aevum/diagnostics/wilson_cycle_lifecycle.py`.
- Added `P82.wilson_cycle_lifecycle_reference` to the tectonics benchmark
  suite.
- P82 creates a scripted Wilson-cycle reference for rift -> spreading ocean ->
  passive/mature ocean -> subduction closure -> arc collision -> suture -> old
  orogen relict.  The benchmark checks persistent basin/lineage identity,
  monotonic phase progression, object-set coverage, gateway causality, and
  suture-to-old-orogen inheritance.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p82_wilson_cycle_lifecycle_reference_20260627/`.
- Current metrics:
  scripted frame count `7`, unique basin/lineage counts `1/1`, required object
  sets observed `8/8`, gateway status count `6`, parent link failures `0`,
  old-orogen relicts `1`, and final scripted basin age `520 Myr`.
- Current generated-world audit:
  ocean basins `7`, Wilson cycles `7`, ocean gateways `11`, rift systems `4`,
  passive margins `2`, closing margins `4`, sutures `1`, and active phase-code
  count `4`.
- Recorded residual:
  current generated Earth-like worlds still lack `tectonics.spreading_centers`.

Next:

- Implement `P83.crust_sediment_province_coupling`.

2026-06-27 - Post-P82 staged research execution plan archived

- Expanded this archive with the post-P82 execution control plan for stages
  0-7.
- The plan keeps climate, ocean circulation, and monsoon redesign paused while
  the terrain/plate system is compared against real-Earth geomorphology.
- The stages now explicitly cover:
  current-state gap inventory, source-corpus provenance, global
  hypsometry/bathymetry/planform envelopes, province-crust-sediment-basement
  coupling, boundary/deep-time lifecycle residuals, mountain/rift/margin
  feature families, drainage/erosion/source-to-sink coupling, and final
  integrated Earth-like promotion audit.
- Each stage records source groups, theory questions, implementation method,
  candidate microbenchmarks, and optimization/acceptance targets so future
  work can proceed one process layer at a time.
- New candidate audit names are reserved for later execution:
  `P90.current_world_morphology_gap_inventory` and
  `P91.integrated_real_earth_morphology_promotion_audit`.

Next:

- Execute Stage 3's immediate implementation entry:
  `P83.crust_sediment_province_coupling`.

2026-06-27 - P83 crust-sediment-province coupling completed

- Added `aevum/diagnostics/crust_sediment_province_coupling.py`.
- Added `P83.crust_sediment_province_coupling` and `run_p83_bench` to the
  tectonics benchmark CLI.
- P83 defines a deterministic synthetic continent reference for province,
  crust, sediment, basement age, elevation, and relief coupling.  It covers
  shield, platform, intracratonic basin, foreland basin, active orogen, old
  orogen, old suture, rift shoulder, rift basin, rift axis, passive-margin
  lowland, continental shelf, and volcanic/LIP plateau classes.
- The reference fixture uses source IDs `CRUST1_0`,
  `GLIM_GLOBAL_LITHOLOGY`, and `NOAA_TOTAL_SEDIMENT_THICKNESS` as provenance
  but stores no raw external rasters or vectors.
- Acceptance now checks crust-thickness ordering, sediment-accommodation
  ordering, elevation and relief ordering, basin/passive-margin lows without
  erasing the parent continent, and the key shield rule: old/stable does not
  mean uniformly high and flat.
- P83 also audits a current 900-cell Earth-like generated world using
  `terrain.continental_landforms`, `crust.thickness_m`,
  `sediment.thickness_m`, and `terrain.elevation_m`.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p83_crust_sediment_province_coupling_20260627/`.
- Reference metrics:
  province classes `13`, parent processes `13`, class edges `23`, missing
  required classes/processes `0/0`, shield elevation/relief `620/380 m`,
  foreland sediment `3600 m`, active-orogen crust thickness `58000 m`, and
  passive-margin lowland elevation `40 m`.
- Current generated-world audit:
  continental landform objects `31`, kind count `8`, basin/lowland mean
  elevation `-222 m` versus platform `1469 m`, basin/lowland mean sediment
  `2095 m` versus platform `703 m`.
- Recorded residual:
  first-class production fields `tectonics.continental_province_id`,
  `tectonics.continental_province_code`, and
  `tectonics.province_parent_process` are still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass -q`
  -> `1 passed in 30.85s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P83 --out out_bench_p83_crust_sediment_province_coupling_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p79_province_reference_graph_extraction_pass tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p81_boundary_process_geometry_reference_pass tests/test_tectonics_bench.py::test_p82_wilson_cycle_lifecycle_reference_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass -q`
  -> `5 passed in 123.70s`.

Next:

- Implement `P84.source_to_sink_sediment_budget`.

2026-06-27 - P84 source-to-sink sediment budget completed

- Added `aevum/diagnostics/source_to_sink_sediment_budget.py`.
- Added `P84.source_to_sink_sediment_budget` and `run_p84_bench` to the
  tectonics benchmark CLI.
- P84 defines a deterministic source-to-sink sediment budget fixture with
  mountain and platform source zones, foreland/passive-margin/shelf/ocean-basin
  sinks, and explicit routing edges.
- The reference checks sediment mass balance, source export, sink deposition,
  routing-edge closure, accommodation utilization, erosion availability, and
  projected land-mask stability before source-to-sink sediment is allowed to
  drive production terrain more strongly.
- P84 also audits a current 900-cell Earth-like generated world and records
  that production source-to-sink objects are still missing.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p84_source_to_sink_sediment_budget_20260627/`.
- Reference metrics:
  zones `6`, routing edges `5`, source/sink volumes `69000/69000 km3`,
  volume-balance fraction `0.0`, max accommodation utilization `0.68`,
  land-mask changes `0`, and routing mismatches `0`.
- Routing volumes:
  mountain export `45000 km3`, platform export `24000 km3`, foreland
  deposition `22500 km3`, passive-margin deposition `18000 km3`, shelf
  deposition `20400 km3`, and ocean-basin deposition `8100 km3`.
- Current generated-world residual:
  `terrain.drainage_basins`, `terrain.sediment_routing_edges`, and
  `terrain.sediment_budget` are still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass -q`
  -> `1 passed in 31.08s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P84 --out out_bench_p84_source_to_sink_sediment_budget_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass -q`
  -> `2 passed in 61.70s`.

Next:

- Implement `P85.drainage_divide_province_alignment`.

2026-06-27 - P85 drainage-divide province alignment completed

- Added `aevum/diagnostics/drainage_divide_province_alignment.py`.
- Added `P85.drainage_divide_province_alignment` and `run_p85_bench` to the
  tectonics benchmark CLI.
- P85 defines a deterministic multi-province continent fixture with shield,
  platform, intracratonic basin, active orogen, foreland basin, rift shoulder,
  rift axis, rift basin, passive-margin lowland, and continental shelf
  provinces.
- The reference fixture checks that drainage divides align with province
  boundaries and highland/divide parent provinces, that flow paths reach
  expected sinks, that paths do not cross divides or other basins, and that
  drainage basins are contiguous rather than checkerboarded.
- P85 also audits a current 900-cell Earth-like generated world and records
  that production drainage objects are still missing.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p85_drainage_divide_province_alignment_20260627/`.
- Reference metrics:
  province classes `10`, drainage basins `3`, divide cells `5`, divide
  fraction `0.0625`, divide alignment `1.0`, highland alignment `1.0`, flow
  paths `6`, flow-to-sink consistency `1.0`, downhill-step fraction `1.0`,
  uphill/divide-crossing/basin-crossing/sink-failure counts `0/0/0/0`, and max
  basin component count `1`.
- Current generated-world residual:
  `terrain.drainage_basins`, `terrain.drainage_divides`,
  `terrain.flow_direction`, and `terrain.flow_accumulation` are still missing.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass -q`
  -> `1 passed in 31.11s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P85 --out out_bench_p85_drainage_divide_province_alignment_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass -q`
  -> `2 passed in 61.83s`.

Next:

- Implement `P86.old_orogen_erosion_decay`.

2026-06-27 - P86 old-orogen erosion decay completed

- Added `aevum/diagnostics/old_orogen_erosion_decay.py`.
- Added `P86.old_orogen_erosion_decay` and `run_p86_bench` to the tectonics
  benchmark CLI.
- P86 defines a deterministic old-orogen decay fixture spanning
  `active_orogen -> post_collision_high_orogen -> decaying_orogen ->
  old_orogen -> subdued_old_orogen`.
- The reference gate checks monotonic relief/elevation/crustal-root decay,
  late sediment-export decline, retained inherited suture/boundary memory, and
  explicit parent-process links.
- P86 also audits a current 900-cell Earth-like generated world and records
  that production old-orogen decay budget fields are still missing.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p86_old_orogen_erosion_decay_20260627/`.
- Reference metrics:
  frame count `5`, relief `1900 -> 420 m`, relief-decay fraction
  `0.7789473684210526`, elevation `3200 -> 680 m`, crustal thickness
  `58000 -> 42000 m`, boundary strength `1.0 -> 0.62`, minimum boundary-trace
  overlap `1.0`, total sediment export `62160 km3`, peak interval sediment
  export `21000 km3`, final interval sediment export `9660 km3`, and parent
  link failures `0`.
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
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P86 --out out_bench_p86_old_orogen_erosion_decay_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass -q`
  -> `2 passed in 61.75s`.

Next:

- Implement `P87.mountain_inventory_expression`.

2026-06-27 - P87 mountain inventory expression completed

- Added `aevum/diagnostics/mountain_inventory_expression.py`.
- Added `P87.mountain_inventory_expression` and `run_p87_bench` to the
  tectonics benchmark CLI.
- P87 defines a small GMBA-style derived reference fixture for mountain ranges
  without bundling raw external mountain vectors.  The fixture covers active
  margin orogens, active collision orogens, collision plateaus, old subdued
  orogens, rift-shoulder ranges, volcanic arc chains, and extensional ranges.
- The reference gate checks range count, mountain class diversity, hierarchy
  levels, parent range links, parent processes, object backing, area caps,
  elongation distribution, and relief envelope.
- P87 also audits a current 900-cell Earth-like generated world and records
  that current mountain expression is object-backed but still incomplete.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p87_mountain_inventory_expression_20260627/`.
- Reference metrics:
  ranges `11`, mountain classes `7`, hierarchy levels `3`, parent processes
  `9`, parent-link failures `0`, threshold-only ranges `0`, total mountain
  area fraction `0.066`, max range area fraction `0.015`, median elongation
  `6.4`, elongated ranges `8`, and max relief `5000 m`.
- Current generated-world audit:
  continental landform objects `31`, mountain candidate objects `11`,
  expressed mountain objects `7`, parented/parent-linked mountain objects
  `11/11`, mountain kind count `2`, total mountain candidate area fraction
  `0.14999228221409389`, expressed mountain area fraction
  `0.11888465637266141`, max mountain object area fraction
  `0.04999122531021068`, parent process/context coverage `1.0/1.0`, median
  elongation `1.8097984605150934`, elongated mountain object count `0`, and
  max relief `4036.162872011885 m`.
- Current generated-world residual:
  `orogen` and `plateau` mountain kinds are still missing; first-class fields
  `terrain.mountain_ranges`, `terrain.mountain_inventory`,
  `terrain.mountain_hierarchy_level`, `tectonics.mountain_belt_id`, and
  `tectonics.mountain_parent_process_id` are still missing; elongated range
  expression is still underdeveloped.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass -q`
  -> `1 passed in 31.39s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P87 --out out_bench_p87_mountain_inventory_expression_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass -q`
  -> `2 passed in 62.00s`.

Next:

- Implement `P88.rift_margin_escarpment_sequence`.

2026-06-27 - P88 rift-margin escarpment sequence completed

- Added `aevum/diagnostics/rift_margin_escarpment_sequence.py`.
- Added `P88.rift_margin_escarpment_sequence` and `run_p88_bench` to the
  tectonics benchmark CLI.
- P88 defines a deterministic rift-to-passive-margin reference transect:
  `stable_platform -> rift_shoulder -> rift_basin -> rift_axis ->
  opposite_rift_shoulder -> passive_margin_escarpment ->
  passive_margin_lowland -> continental_shelf -> continental_slope ->
  continental_rise -> abyssal_plain`.
- The reference gate checks required adjacency, rift shoulder/basin relief
  ordering, escarpment expression, passive-margin lowland/shelf coupling,
  shelf-slope-rise-abyss depth ordering, sediment ordering, and parent-process
  coverage.
- P88 also audits a current 900-cell Earth-like generated world and records
  that rift/passive-margin/shelf coupling is present, but first-class
  rift-margin sequence objects are still missing.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p88_rift_margin_escarpment_sequence_20260627/`.
- Reference metrics:
  zones `11`, classes `10`, adjacency edges `10`, missing required edges `0`,
  parent process failures `0`, rift shoulder/basin elevations `980/-80 m`,
  escarpment relief `780 m`, passive-margin lowland elevation `55 m`, shelf /
  slope / rise / abyss depths `120/1700/3100/4300 m`, shelf sediment `3200 m`,
  passive-margin lowland sediment `1850 m`, and rift-basin sediment `2800 m`.
- Current generated-world audit:
  continental landform objects `31`, margin landform objects `16`, rift-basin
  objects `7`, passive-margin lowland objects `2`, passive-margin wedge objects
  `2`, delta/fan objects `2`, rift-basin area fraction
  `0.07669787118898672`, passive-margin lowland object/field area fractions
  `0.004445749597838481/0.007779671667038456`, passive-margin wedge area
  fraction `0.22223541278731013`, shelf/slope/rise/abyss ocean fractions
  `0.09085636887296934/0.10224581537585933/0.0876907663053517/0.3083919118087776`,
  lowland-near-shelf fraction at two passes `1.0`, lowland-near-wedge fraction
  at two passes `1.0`, rift-near-passive-margin fraction at five passes
  `0.985508126394939`, shelf depth p75 `122.88691954044891 m`, abyss depth p50
  `3509.9321625103657 m`, and parented rift/lowland/wedge objects `7/2/2`.
- Current generated-world residual:
  `terrain.rift_shoulders`, `terrain.escarpments`,
  `terrain.rift_margin_sequence_id`, `terrain.rift_margin_stage`, and
  `tectonics.rift_margin_lineage_id` are still missing; passive-margin lowland
  objects remain tiny.
- Verification:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass -q`
  -> `1 passed in 31.11s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P88 --out out_bench_p88_rift_margin_escarpment_sequence_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p87_mountain_inventory_expression_pass tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass -q`
  -> `2 passed in 61.40s`.

Next:

- Implement `P89.plateau_area_cap_and_decay`.

2026-06-27 - Staged evidence collection and benchmark enrichment archive added

- Added the `Staged Evidence Collection and Benchmark Enrichment Archive`
  section.
- The new archive defines a reusable evidence-packet contract so future source
  collection produces source metadata, extraction notes, theory claims,
  derived metrics, reference fixtures, generated-world audits, optimization
  targets, residual policy, and asset-review requirements.
- Added staged research batches:
  - `R0.current_generated_world_forensic_baseline`;
  - `R1.global_topography_bathymetry_planform_evidence`;
  - `R2.continental_province_basement_crust_sediment_evidence`;
  - `R3.plate_boundary_wilson_cycle_deep_time_evidence`;
  - `R4.mountains_plateaus_rifts_margins_special_landform_evidence`;
  - `R5.drainage_erosion_source_to_sink_evidence`;
  - `R6.real_earth_case_study_evidence_packets`;
  - `R7.integrated_calibration_and_promotion_matrix`.
- The archive explicitly ties each batch to theory basis, implementation
  method, microbenchmark enrichment, and optimization targets.
- The current execution order remains unchanged: finish `P89` first, then use
  the expanded archive to define `P90.current_world_morphology_gap_inventory`
  and `P91.integrated_real_earth_morphology_promotion_audit`.

Next:

- Implement `P89.plateau_area_cap_and_decay`, then use the new evidence-packet
  contract to scope `P90`.

2026-06-27 - P89 plateau area cap and decay completed

- Added `aevum/diagnostics/plateau_area_cap_and_decay.py`.
- Added `P89.plateau_area_cap_and_decay` and `run_p89_bench` to the tectonics
  benchmark CLI.
- P89 defines a deterministic two-variant plateau lifecycle reference:
  collision plateaus from incipient collision through mature/post-peak/
  dissected stages, and volcanic/LIP plateaus from plume swell through LIP
  emplacement/post-LIP/eroded-surface stages.
- The reference gate requires parent processes and parent object kinds, finite
  plateau area caps, background-interior dominance, collision plateau decay,
  volcanic plateau decay, and plateau elevation above surrounding platforms.
- P89 also audits a current 900-cell Earth-like generated world and records
  that plateau area overpaint is currently absent, but production plateau
  expression is still missing.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p89_plateau_area_cap_and_decay_20260627/`.
- Reference metrics:
  variants/frames/stages `2/8/8`, max collision plateau area fraction
  `0.024`, max volcanic plateau area fraction `0.016`, combined peak plateau
  area fraction `0.04`, collision elevation decay `2450 m`, collision area
  decay `0.013`, volcanic elevation decay `1050 m`, volcanic area decay
  `0.01`, and parent process/object failures `0/0`.
- Current generated-world audit:
  plateau objects `0`, plateau detail cells `0`, LIP objects `8`, first-class
  plateau fields missing, expected plateau kinds `plateau` and
  `volcanic_lip_plateau` missing, and high interior without plateau support
  fraction `0.08958631689930298` of continental land.
- Current generated-world residual:
  `terrain.plateau_inventory`, `terrain.plateau_age_myr`,
  `terrain.plateau_decay_stage`, `terrain.plateau_parent_process_id`, and
  `tectonics.plateau_lineage_id` are still missing; plateau and volcanic/LIP
  plateau expression remain production gaps.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/plateau_area_cap_and_decay.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p89_plateau_area_cap_and_decay_pass -q`
  -> `1 passed in 31.23s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P89 --out out_bench_p89_plateau_area_cap_and_decay_20260627`
  -> `status: pass`.
  P88-P89 chain:
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p88_rift_margin_escarpment_sequence_pass tests/test_tectonics_bench.py::test_p89_plateau_area_cap_and_decay_pass -q`
  -> `2 passed in 63.47s`.

Next:

- Implement `P90.current_world_morphology_gap_inventory` from the R0 evidence
  packet contract, using P76-P89 outputs to group current defects by owner
  layer.

2026-06-27 - P90 current-world morphology gap inventory completed

- Added `aevum/diagnostics/current_world_morphology_gap_inventory.py`.
- Added `P90.current_world_morphology_gap_inventory` and `run_p90_bench` to the
  tectonics benchmark CLI.
- P90 builds one 900-cell Earth-like generated world, reuses P78-P89 current
  audits plus compiler consistency, and groups current defects by owner layer
  and category.
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
- Key residual groups:
  planform calibration, missing `volcanic_lip_plateau`, missing transform and
  `tectonics.spreading_centers`, missing production province fields,
  drainage/source-to-sink fields, old-orogen decay fields, mountain/plateau
  inventory fields, rift-margin sequence fields, and missing P91 review PNG
  set.
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

2026-06-27 - P91 integrated real-Earth morphology promotion audit completed

- Added `aevum/diagnostics/integrated_real_earth_morphology_promotion_audit.py`.
- Added `P91.integrated_real_earth_morphology_promotion_audit` and
  `run_p91_bench` to the tectonics benchmark CLI.
- P91 reads archived P76-P90 benchmark summaries, verifies the archived P69
  8000-cell high-resolution PNG/contact-sheet evidence, generates fresh
  900-cell and 2500-cell Earth-like CI world assets, attaches those assets to
  the current-world gap inventory, and records a promotion decision.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p91_integrated_real_earth_morphology_promotion_audit_20260627/`.
- Generated P91 review assets:
  `p91_ci_world_contact_sheet.png`, `ci_world_assets/earthlike_900cells/`,
  and `ci_world_assets/earthlike_2500cells/`.
- Current P91 metrics:
  P76-P90 stage pass count `15/15`, archived high-resolution required PNGs
  `24/24`, high-resolution contact sheets `1`, CI required PNGs `16/16`, CI
  compilers passed `2/2`, CI inventories ready `2/2`, root P90 non-asset gaps
  `50`, owner layers `7`, and unassigned/generic blockers `0/0`.
- Promotion decision:
  `audit_completed=True`, `promotion_ready=False`, and
  `promotion_decision_recorded=True`.
- Named promotion blockers:
  `p69_earthlike_reference_needs_calibration`,
  `p90_current_world_residuals_unresolved`,
  `bathymetry_margin_residuals_unresolved`,
  `boundary_lifecycle_residuals_unresolved`,
  `crust_sediment_residuals_unresolved`,
  `drainage_erosion_residuals_unresolved`,
  `landform_expression_residuals_unresolved`,
  `planform_residuals_unresolved`, and
  `province_graph_residuals_unresolved`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/integrated_real_earth_morphology_promotion_audit.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_integrated_real_earth_morphology_promotion_audit_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p91_integrated_real_earth_morphology_promotion_audit_pass -q`
  -> `1 passed in 130.21s`.

Next:

- Define and implement `P92.production_residual_owner_repair_plan`, starting
  from the P91 blocker matrix rather than reopening climate/ocean-current work.

2026-06-27 - P92 production residual owner repair plan completed

- Added `aevum/diagnostics/production_residual_owner_repair_plan.py`.
- Added `P92.production_residual_owner_repair_plan` and `run_p92_bench` to the
  tectonics benchmark CLI.
- P92 consumes the latest archived P91 audit rather than rerunning the
  expensive 900/2500-cell render pass.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p92_production_residual_owner_repair_plan_20260627/`.
- Current P92 metrics:
  P91 blockers `9/9` assigned, owner layers `7/7` assigned, residual items
  `32/32` assigned, repair packets `8`, dependency order valid, climate/ocean/
  monsoon targets `0`, and final P91 reaudit defined.
- Repair packet sequence:
  `P92.1_planform_and_reference_calibration`,
  `P92.2_production_province_graph_fields`,
  `P92.3_boundary_lifecycle_objects`,
  `P92.4_crust_sediment_interior_relief_coupling`,
  `P92.5_drainage_source_to_sink_fields`,
  `P92.6_landform_inventory_lifecycle`,
  `P92.7_bathymetry_margin_sequence`, and
  `P92.8_integrated_reaudit_and_promotion_gate`.
- Next implementation packet:
  `P92.1_planform_and_reference_calibration`, with candidate microbenchmarks
  `P93.planform_reference_calibration` and
  `P93.generated_component_ribbon_envelope`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/production_residual_owner_repair_plan.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P92 --out out_bench_p92_production_residual_owner_repair_plan_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p92_production_residual_owner_repair_plan_pass -q`
  -> `1 passed in 0.87s`.

Next:

- Implement `P93.planform_reference_calibration`, starting with the
  planform/reference blocker before province graph and downstream terrain
  repairs.

2026-06-27 - P93 planform reference calibration completed

- Added `aevum/diagnostics/planform_reference_calibration.py`.
- Added `P93.planform_reference_calibration`,
  `P93.generated_component_ribbon_envelope`, and `run_p93_bench` to the
  tectonics benchmark CLI.
- P93 consumes archived evidence only: P69 high-resolution physical ensemble,
  P78 generated hypsometry envelope, P90 current-world gap inventory, P91
  integrated promotion audit, and P92 production residual owner repair plan.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p93_planform_reference_calibration_20260627/`.
- Current P93 metrics:
  P90 planform gaps `9`, calibration targets `5`, covered planform metrics
  `5`, cross-owner deferred targets `1`, P69 reference members `3`, P69
  Earth-like out-of-envelope count `6`, P78 current out-of-envelope count `5`,
  unresolved primary planform metrics `5`, and next packet
  `P92.2_production_province_graph_fields`.
- Planform targets archived:
  exposed land fraction too low, major land component count too low,
  exposed-land ribbon fraction too high, largest-landmass coastline complexity
  too high, and high-resolution Earth-like largest-component share too high.
- Cross-owner target archived:
  trench fraction is too high, but its production owner is
  `P92.7_bathymetry_margin_sequence`, not the P94 province-graph packet.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/planform_reference_calibration.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P93 --out out_bench_p93_planform_reference_calibration_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p93_planform_reference_calibration_pass -q`
  -> `1 passed in 0.79s`.

Next:

- Implement `P94.production_province_graph_fields`, the production benchmark for
  `P92.2_production_province_graph_fields`.

2026-06-27 - P94 production province graph fields completed

- Added the production continental province graph layer:
  `terrain.continental_province_id`, `terrain.continental_province_code`,
  `tectonics.continental_province_id`,
  `tectonics.continental_province_code`,
  `tectonics.province_parent_process`, and
  `tectonics.continental_provinces`.
- Added P94 diagnostics for 900/2500-cell generated worlds:
  production field coverage, object consistency, parent-process coverage,
  multi-province major continents, required class coverage, LIP/rift adjacency,
  and checkerboard suppression.
- Updated P80/P83 so the generated-world reference comparison and
  crust-sediment province audit consume the production province graph instead
  of recording it as a pending residual.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p94_production_province_graph_fields_20260627/`.
- Current P94 metrics:
  ready worlds `2/2`, minimum province object count `43`, minimum field id
  count `43`, minimum province class count `9`, minimum parent-process count
  `10`, id/code/parent-process coverage `1.0`, max missing field/class/edge
  count `0`, max object mismatch count `0`, max disconnected province id count
  `0`, max tiny-province area fraction `0.01664740312578796`, minimum
  `rift_system|volcanic_lip_plateau` edge count `2`, and P80 900-cell residual
  class/edge counts `0`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/terrain.py aevum/diagnostics/production_province_graph.py aevum/diagnostics/tectonics_bench.py aevum/diagnostics/generated_province_reference.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P94 --out out_bench_p94_production_province_graph_fields_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p94_production_province_graph_fields_pass -q`
  -> `3 passed in 179.99s`.

Next:

- Implement `P95.boundary_lifecycle_objects`, the production benchmark for
  `P92.3_boundary_lifecycle_objects`.

2026-06-27 - P95 boundary lifecycle objects completed

- Added production transform boundary support:
  ridge-offset transform process geometry is generated when neutral shear is
  under-resolved after boundary thinning.
- Added sparse ridge/transform boundary-object aggregation so lifecycle
  objects survive spherical grid thinning instead of being dropped as single
  cells.
- Updated P81/P82 diagnostics and tests so transform process geometry and
  `tectonics.spreading_centers` are required current production outputs rather
  than expected residuals.
- Added P95 diagnostics for 900/2500-cell generated worlds:
  transform process coverage, ridge/transform boundary objects,
  spreading-center lifecycle objects, current boundary/Wilson readiness,
  phase diversity, gateway counts, and parent-link integrity.
- Current output:
  `/Users/rayw/Projects/aevum/out_bench_p95_boundary_lifecycle_objects_20260627/`.
- Current P95 metrics:
  ready worlds `2/2`, minimum transform cell count `15`, minimum ridge cell
  count `23`, minimum transform boundary object count `3`, minimum ridge
  boundary object count `4`, minimum spreading-center count `4`, minimum
  transform-near-ridge fraction `0.9655172413793104`, max transform length
  fraction `0.08064103031499571`, missing boundary process type count `0`,
  missing Wilson object-set count `0`, and parent link failures `0`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/tectonics.py aevum/diagnostics/boundary_process_geometry.py aevum/diagnostics/wilson_cycle_lifecycle.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P95 --out out_bench_p95_boundary_lifecycle_objects_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p81_boundary_process_geometry_reference_pass tests/test_tectonics_bench.py::test_p82_wilson_cycle_lifecycle_reference_pass tests/test_tectonics_bench.py::test_p95_boundary_lifecycle_objects_pass -q`
  -> `3 passed in 177.16s`.

Next:

- Implement `P96.high_flat_interior_owner_reduction` and
  `P96.province_crust_sediment_surface_ordering` for
  `P92.4_crust_sediment_interior_relief_coupling`.

2026-06-27 - P96 crust/sediment interior relief coupling completed

- Added production surface-ordering logic in terrain:
  crust thickness/support, stable platform state, passive margins, foreland
  basins, rift potential, sediment accommodation, and deformation context now
  adjust exposed continental elevation before the final generated-world audits.
- Updated generated-world crust/sediment coupling:
  when `tectonics.continental_provinces` exists, current-world metrics are
  grouped by production province classes instead of legacy
  `terrain.continental_landforms`.
- Updated generated province-reference comparison:
  production province fields are now authoritative for generated class
  coverage, while legacy overlays remain counted as overlays rather than
  overwriting production classes.
- Updated production province graph coverage:
  deterministic anchors now preserve platform, shield, passive-margin lowland,
  intracratonic basin, and foreland-basin representation across major
  generated continents; deliberate anchors are excluded from the
  checkerboard-noise fraction.
- Added P96 diagnostics for 900/2500-cell generated worlds:
  high-flat interior owner reduction, basin/lowland preservation,
  production-graph readiness, production coupling aggregation, P80/P94
  residual clearance, and platform/basin/orogen elevation/sediment ordering.
- Current outputs:
  `/Users/rayw/Projects/aevum/out_bench_p96_crust_sediment_surface_ordering_20260627_d/`.
  `/Users/rayw/Projects/aevum/out_bench_p94_after_p96_20260627_b/`.
- Current P96 metrics:
  ready worlds `2/2`, maximum direct high-flat interior fraction
  `0.06246298636795882`, owner-aware P96 high-flat fraction after ordering
  `0.0`, minimum basin/lowland fraction `0.37219087217151403`, minimum
  `<1000m` lowland fraction `0.5185542116145031`, minimum production province
  class count `9`, minimum parent-process count `10`, P80/P94 missing class and
  edge counts `0`, basin-minus-platform elevation at most
  `-686.6328789803508m`, basin-minus-platform sediment at least
  `941.43280342115m`, and orogen-minus-basin elevation at least
  `554.1399306870508m`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/terrain.py aevum/diagnostics/crust_sediment_province_coupling.py aevum/diagnostics/generated_province_reference.py aevum/diagnostics/production_province_graph.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
  -> passed.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P96 --out out_bench_p96_crust_sediment_surface_ordering_20260627_d`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P94 --out out_bench_p94_after_p96_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p80_generated_province_graph_reference_comparison_pass tests/test_tectonics_bench.py::test_p83_crust_sediment_province_coupling_pass tests/test_tectonics_bench.py::test_p94_production_province_graph_fields_pass tests/test_tectonics_bench.py::test_p96_crust_sediment_surface_ordering_pass -q`
  -> `4 passed in 306.66s`.

Next:

- Implement `P97.production_drainage_source_to_sink_fields` for
  `P92.5_drainage_source_to_sink_fields`.

2026-06-27 - P97 drainage/source-to-sink production fields completed

- Added production drainage fields:
  `terrain.drainage_basins`, `terrain.drainage_divides`,
  `terrain.flow_direction`, `terrain.flow_accumulation`, and
  `terrain.drainage_surface_m`.
- Added production source-to-sink sediment fields and objects:
  `terrain.sediment_source_m`, `terrain.sediment_sink_m`,
  `terrain.sediment_budget_balance`, basin objects, divide objects, sediment
  routing edges, and a sediment budget object.
- Added production old-orogen lifecycle/erosion fields:
  `terrain.old_orogen_decay_stage`,
  `terrain.orogen_erosion_budget`,
  `terrain.orogen_boundary_memory`, and
  `terrain.orogen_sediment_export`.
- Updated P84/P85/P86 so the current generated-world audits require these
  production outputs instead of treating them as unresolved residuals.
- Added P97 diagnostics for 900/2500-cell generated worlds and updated P90/P91
  reaudit logic so current same-run summaries replace stale archived stage
  summaries.
- Current outputs:
  `/Users/rayw/Projects/aevum/out_bench_p97_drainage_source_to_sink_fields_20260627/`.
  `/Users/rayw/Projects/aevum/out_bench_p90_after_p97_probe_20260627_b/`.
  `/Users/rayw/Projects/aevum/out_bench_p91_after_p97_probe_20260627_b/`.
- Current P97 metrics:
  ready worlds `2/2`, minimum drainage basin object count `18`, minimum routing
  edge count `39`, minimum routing source/sink kind counts `3/5`, maximum
  source-to-sink balance fraction `1.8436625312242433e-16`, maximum drainage
  divide fraction of land `0.5095258030376043`, divide alignment `1.0`,
  major basin component failures `0`, flow-to-sink consistency `1.0`,
  downhill path fraction `1.0`, old-orogen decay stage count at least `2`,
  minimum boundary memory `0.6682772168847143`, minimum erosion budget
  `38.79285174106078m`, and minimum sediment export
  `2239589.12726769km3`.
- Verification:
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P97 --out out_bench_p97_drainage_source_to_sink_fields_20260627`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p84_source_to_sink_sediment_budget_pass tests/test_tectonics_bench.py::test_p85_drainage_divide_province_alignment_pass tests/test_tectonics_bench.py::test_p86_old_orogen_erosion_decay_pass tests/test_tectonics_bench.py::test_p97_drainage_source_to_sink_fields_pass -q`
  -> `4 passed in 228.06s`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P90 --out out_bench_p90_after_p97_probe_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_after_p97_probe_20260627_b`
  -> `status: pass`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p90_current_world_morphology_gap_inventory_pass -q`
  -> `1 passed in 34.90s`.
  `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p91_integrated_real_earth_morphology_promotion_audit_pass -q`
  -> `1 passed in 155.84s`.

Next:

- Implement `P98.production_landform_inventory_lifecycle_fields` for
  `P92.6_landform_inventory_lifecycle`.

2026-06-27 - P98 landform inventory/lifecycle production fields completed

- Added production mountain inventory fields:
  `terrain.mountain_ranges`, `terrain.mountain_inventory`,
  `terrain.mountain_hierarchy_level`, `tectonics.mountain_belt_id`, and
  `tectonics.mountain_parent_process_id`.
- Added production plateau inventory/lifecycle fields:
  `terrain.plateau_inventory`, `terrain.plateau_age_myr`,
  `terrain.plateau_decay_stage`, `terrain.plateau_parent_process_id`, and
  `tectonics.plateau_lineage_id`.
- Added object-backed `terrain.mountain_ranges` and
  `terrain.plateau_inventory` production layers and extended
  `terrain.continental_landforms` with process-parented P98 orogen/plateau
  objects.
- Updated P87/P89 so current generated-world audits require these production
  fields and no longer report mountain/plateau inventory gaps as current
  residuals.
- Updated P90/P91 reaudit logic so same-run P87/P89 summaries replace stale
  archived mountain/plateau residual state.
- Current outputs:
  `/Users/rayw/Projects/aevum/out_bench_p98_landform_inventory_lifecycle_20260627_d/`.
  `/Users/rayw/Projects/aevum/out_bench_p90_after_p98_probe_20260627_b/`.
  `/Users/rayw/Projects/aevum/out_bench_p91_after_p98_probe_20260627_b/`.
- Current P98 metrics:
  ready worlds `2/2`, mountain missing field/kind counts `0/0`, minimum
  production mountain range object count `6`, minimum mountain field id count
  `6`, minimum mountain inventory class count `4`, minimum maximum mountain
  elongation ratio `1.748363423468797`, maximum mountain area fraction
  `0.1244333061118279`, plateau missing item/kind counts `0/0`, minimum
  plateau inventory cell count `19`, minimum volcanic/LIP plateau cell count
  `19`, maximum plateau area fraction `0.021075501103753724`, and maximum
  high-interior-without-plateau fraction `0.07413285329377038`.
- Reaudit state:
  P90 after-P98 records `25` gaps, `13` current residual items, `4` owner
  layers, and no missing mountain/plateau production fields.
  P91 after-P98 records `15/15` stage summaries passing, CI assets `16/16`,
  CI inventories `2/2`, root non-asset gaps `17`, root owner layers `3`, and
  promotion blockers `5`: P69 reference calibration, unresolved P90 residuals,
  bathymetry/margin, crust/sediment, and planform.
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

2026-06-27 - P99 bathymetry/margin sequence production fields completed

- Production work:
  - Added first-class rift-margin sequence fields:
    `terrain.rift_margin_sequence_id`, `terrain.rift_margin_stage`, and
    `tectonics.rift_margin_lineage_id`.
  - Added first-class rift-shoulder and escarpment expression:
    `terrain.rift_shoulders`, `terrain.escarpments`, and corresponding
    production objects.
  - Added `terrain.rift_margin_sequences` objects to preserve the linked
    shoulder -> rift basin -> escarpment/passive lowland -> shelf/slope/rise/
    abyss lifecycle as an auditable margin sequence.
- Current outputs:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P99 --out out_bench_p99_bathymetry_margin_sequence_20260627_c`
    -> `status: pass`.
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P91 --out out_bench_p91_after_p99_probe_20260627_b`
    -> `status: pass`.
- Current P99 metrics:
  - `sequence_ready_world_count: 2`, `max_missing_sequence_item_count: 0`,
    `min_sequence_id_count: 1`, `min_lineage_id_count: 1`,
    `min_stage_count: 8`.
  - `min_rift_shoulder_cell_count: 143`,
    `min_escarpment_cell_count: 3`,
    `min_sequence_object_count: 1`,
    `min_rift_shoulder_object_count: 1`,
    `min_escarpment_object_count: 1`.
  - `ordered_world_count: 2`, `min_shelf_stage_cell_count: 62`,
    `min_slope_stage_cell_count: 65`, `min_rise_stage_cell_count: 100`,
    `min_abyss_stage_cell_count: 119`,
    `min_shelf_to_abyss_depth_delta_m: 4083.49`.
- Current after-P99 promotion audit:
  - P90 after-P99: `19` gaps, `3` owner layers, `3` categories,
    `8` current residual items.
  - P91 after-P99: `11` non-asset root gaps, `2` root owner layers,
    `4` blockers, promotion remains default-off.
  - Remaining blockers are `p69_earthlike_reference_needs_calibration`,
    `p90_current_world_residuals_unresolved`,
    `crust_sediment_residuals_unresolved`, and
    `planform_residuals_unresolved`.
- Verification:
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py -k 'p88_rift_margin_escarpment_sequence or p90_current_world_morphology_gap_inventory or p91_integrated_real_earth_morphology_promotion_audit or p99_bathymetry_margin_sequence'`
    -> `4 passed, 88 deselected, 2 warnings in 359.53s`.

Next at that point:

- Implement `P100.integrated_reaudit_and_promotion_gate` for
  `P92.8_integrated_reaudit_and_promotion_gate`, then decide whether the next
  production packet should prioritize planform repair, crust/sediment
  calibration, or a release-promotion gate.

2026-06-27 - P100 integrated reaudit and promotion gate completed

- Diagnostic work:
  - Added `P100.integrated_owner_repair_reaudit` for the P93-P99 repair-suite
    matrix and fresh after-P99 P91 evidence.
  - Added `P100.default_promotion_decision_gate` for the default Earth-like
    promotion decision.
  - Added `test_p100_integrated_reaudit_and_promotion_gate_pass`.
- Current output:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P100 --out out_bench_p100_integrated_reaudit_and_promotion_gate_20260627_c`
    -> `status: pass`.
  - Nested P91 evidence is written under
    `out_bench_p100_integrated_reaudit_and_promotion_gate_20260627_c/p91_reaudit/`.
- Current P100 metrics:
  - `repair_suite_count: 7`, `repair_suite_pass_count: 7`,
    `missing_repair_suite_count: 0`, `failing_repair_suite_count: 0`.
  - `p91_after_p99_status: pass`, `p91_stage_suite_pass_count: 15`,
    `p91_ci_asset_set_complete_count: 2`,
    `p91_ci_compiler_passed_count: 2`.
  - `root_p90_non_asset_gap_count: 11`,
    `root_p90_owner_layer_count: 2`,
    `root_p90_residual_item_count: 0`.
  - `cleared_root_owner_blocker_count: 5`,
    `remaining_owner_blocker_count: 2`,
    `promotion_blocker_count: 4`.
- Current promotion decision:
  - Default promotion is still blocked.
  - Cleared root owner blockers are bathymetry/margin, boundary lifecycle,
    drainage/erosion, landform expression, and province graph.
  - Remaining blockers are `p69_earthlike_reference_needs_calibration`,
    `p90_current_world_residuals_unresolved`,
    `crust_sediment_residuals_unresolved`, and
    `planform_residuals_unresolved`.
  - Next action is `P101.planform_crust_sediment_residual_repair`.
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p100_integrated_reaudit_and_promotion_gate_pass -q`
    -> `1 passed, 2 warnings in 153.22s`.

Next:

- Define and implement `P101.planform_crust_sediment_residual_repair`.

2026-06-27 - P101 Phase 0 current residual attribution completed

- Diagnostic work:
  - Added `P101.planform_residual_baseline`.
  - Added `P101.crust_sediment_high_flat_repair`.
  - Added `test_p101_current_residual_attribution_phase0_pass`.
  - P101 is intentionally diagnostic-only for Phase 0: it reproduces and
    attributes current blockers, and records repair targets without relaxing the
    Earth-like envelopes.
- Current output:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P101 --out out_bench_p101_phase0_current_residual_attribution_20260627`
    -> `status: pass`.
- Current P101 root metrics:
  - Root 900-cell gap state: `19` gaps, `11` non-asset gaps, `8` asset-review
    gaps.
  - Root owner split: `10` planform gaps and `1` crust/sediment gap.
  - Attribution: `11/11` expected non-asset gap ids present, `11/11` exact
    code-owner attributions, `0` fallback attributions.
  - Direct residual values: land fraction `0.23110909068830526`, major component
    count `2`, high-flat interior share `0.06246298636795882`, and basin/lowland
    share `0.43297479705190384`.
  - 2500-cell cross-check records `18` gaps, `10` non-asset gaps, owner split
    `7` planform, `1` crust/sediment, and `2` landform-expression gaps.  This
    confirms that Phase 1-3 should keep multi-resolution checks because some
    residual owners vary with scale.
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p101_current_residual_attribution_phase0_pass -q`
    -> `1 passed in 134.07s`.

Next:

- Start Phase 1 reference evidence expansion with a P102/P101-successor packet
  that enriches P76/P77/P79 source fixtures for planform, island/ribbon,
  coastline, province, crust/sediment, and multi-resolution calibration.
- Then start Phase 2 production planform repair against the P101 baseline,
  without changing the current reference envelopes to make the gate pass.

2026-06-27 - P102 Phase 1 reference evidence packets completed

- Diagnostic work:
  - Added `aevum.diagnostics.reference_evidence_packets`.
  - Added `P102.reference_evidence_packet_matrix`.
  - Added `test_p102_reference_evidence_packet_matrix_pass`.
- Current output:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P102 --out out_bench_p102_reference_evidence_packets_20260627`
    -> `status: pass`.
- Current P102 metrics:
  - Evidence packets: `6`.
  - Source ids covered: `17`.
  - Derived metric keys covered: `38`.
  - Covered metric groups: drainage, hypsometry, ocean bathymetry, planform,
    process parentage, and province architecture.
  - Residual owner coverage: bathymetry/margin, boundary lifecycle,
    compiler/render, crust/sediment, drainage/erosion, landform expression,
    planform, and province graph.
  - Upstream fixture status: P76 source ledger ready, P77 hypsometry fixture
    ready, and P79 province reference graph ready.
  - Packet ids:
    `R1_global_hypsometry_planform`,
    `R2_province_crust_sediment_basement`,
    `R3_boundary_wilson_deeptime`,
    `R4_drainage_erosion_source_to_sink`,
    `R5_landform_margins_mountains_plateaus`, and
    `R6_case_study_feature_catalog`.
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/reference_evidence_packets.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/python -m pytest tests/test_tectonics_bench.py::test_p102_reference_evidence_packet_matrix_pass -q`
    -> `1 passed in 0.98s`.

Next:

- Start Phase 2 production planform repair.  The first implementation target
  should be a P103-style planform mechanism repair gate that consumes the P101
  residual baseline and P102 evidence packets before changing land exposure,
  continent component balance, ribbon pruning, or coastline smoothing.

2026-06-27 - P103+ source-corpus enrichment and repair planning archived

- Planning changes:
  - Added `P103+ Source-Corpus Enrichment and Repair Planning Archive`.
  - The archive defines source tiers for global relief/planform, bathymetry,
    crust/sediment/lithology, plate boundaries, deep-time reconstructions,
    active deformation, drainage, erosion, mountain inventories, and
    continent/ocean case-study packages.
  - It defines staged source collection phases S0-S7: ledger hardening, global
    relief/planform extraction, crust/sediment/province extraction, boundary and
    Wilson-cycle extraction, surface-process/mountain extraction, case-study
    packet expansion, production repair mapping, and final promotion/residual
    audit.
  - It reserves P103-P109 benchmark names for planform repair,
    crust/sediment/interior elevation repair, province-boundary expression,
    boundary/Wilson-cycle consistency, drainage/mountain/source-to-sink
    coupling, real-Earth case-study family coverage, and multi-resolution asset
    promotion audit.
  - It records theory-note requirements and optimization discipline so later
    production patches cannot clear gates by random major-event selection,
    visual-only raster flips, or threshold relaxation without source evidence.
- Scope:
  - This is planning and archival work only.  It does not change terrain,
    tectonics, diagnostics code, or test thresholds.
  - Climate, ocean-current, and monsoon redesign remain paused.

Next:

- Start `P103.planform_mechanism_repair` by converting the P103+ planform
  section into an executable benchmark gate, then make production changes only
  against the reproduced P101 residual baseline and P102 evidence packets.

2026-06-27 - P104A continental mosaic expression first slice completed

- Production change:
  - Added a P104A object-expression pass in `TerrainModule` after production
    continental province graph construction and before landform inventory
    export.
  - The pass maps process-backed `terrain.continental_province_code` classes
    back into `terrain.continental_detail` and `terrain.province`, then applies
    modest province-specific relief:
    shields/platforms remain low-to-moderate, basins and passive/rift lowlands
    are lowered but kept exposed where already land, old/active orogens and
    LIP plateaus gain controlled relief.
  - Added a mature-platform split for overlarge platform interiors.  The split
    uses inherited crust stability/thickness and sediment/low-relief signals to
    create shield-like and sag-basin-like subareas.  It does not introduce
    random terrain texture.
- Test:
  - Added `test_p104a_continental_mosaic_expression_fixture`.
  - Verification:
    `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/terrain.py tests/test_tectonics_bench.py`
    -> pass.
  - Verification:
    `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104a_continental_mosaic_expression_fixture -q`
    -> `1 passed in 0.57s`.
- Generated-world validation:
  - 900-cell seeds `42` and `31415`: no failing major-component province
    diversity items after P104A, with exposed continental detail class counts
    `5` and `6`.
  - 2500-cell seed `31415` after platform split:
    land fraction `0.2704`, major components `2`,
    min province classes per major `5`,
    max largest internal province fraction `0.6830`,
    failing major component count `0`.
  - 2500-cell seed `27182` before platform split already had failing major
    component count `0`, min province classes per major `3`, max largest
    internal province fraction `0.6351`.
  - 2500-cell seed `16180` after platform split:
    land fraction `0.2692`, major components `2`,
    min province classes per major `3`,
    max largest internal province fraction `0.7080`,
    failing major component count `0`.
  - Updated rendered assets for the two previously failing 2500-cell worlds:
    `out_p104a_continental_mosaic_v2_20260627/earthlike_seed31415/` and
    `out_p104a_continental_mosaic_v2_20260627/earthlike_seed16180/`.
- Residuals:
  - P104A first slice improves continental interior variety but does not close
    the wider plate/terrain goal.  The same validation runs still report
    planform/ribbon and coastline-complexity warnings.
  - `generated_world_province_diversity_summary` can still report
    `generated_world_gate_incomplete` because legacy acceptance also expects
    passive-margin lowland object counts; the per-major-component diversity
    failure count is now `0` for the checked P104A worlds.

Next:

- Add an explicit P104A generated-world bench summary that separates
  continental-interior mosaic acceptance from legacy passive-margin object
  acceptance.
- Continue P103/P104A iteration on remaining planform/ribbon/coastline issues
  and then proceed to P104B crust/sediment elevation retuning.

2026-06-27 - P104A generated-world mosaic gate completed

- Diagnostic work:
  - Added suite `P104A` to `aevum.diagnostics.tectonics_bench`.
  - Added `P104A.generated_world_continental_mosaic_gate` for multi-seed
    generated-world acceptance of major-continent internal province diversity.
  - Added `P104A.legacy_gate_separation` so the plan distinguishes
    continental-interior mosaic readiness from older passive-margin lowland
    object expectations.
  - Added `test_p104a_generated_world_continental_mosaic_gate_pass`.
- Current output:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P104A --out out_bench_p104a_continental_mosaic_gate_20260627`
    -> `status: pass`, deterministic `True`.
  - Summary and CSV:
    `out_bench_p104a_continental_mosaic_gate_20260627/tectonics_bench_summary.json`
    and
    `out_bench_p104a_continental_mosaic_gate_20260627/p104a_microbenchmarks.csv`.
- Acceptance evidence:
  - Seeds checked: `42`, `31415`, `16180`; resolution `900` cells.
  - Per-major-component internal mosaic pass count: `3/3`.
  - Minimum province classes per major component: `4`.
  - Maximum largest internal province fraction: `0.6430`.
  - Minimum exposed-continental detail class count from P104A telemetry: `6`.
  - Minimum exposed-continental province class count from P104A telemetry: `9`.
  - Minimum basin/lowland share per major component: `0.6369`.
  - Maximum active highland/plateau fraction: `0.1714`.
  - Maximum unparented highland fraction: `0.0`.
  - Mean absolute P104A surface delta remains moderate, max `256.27 m`.
- Interpretation:
  - P104A now proves that current generated worlds are no longer failing the
    narrow "single broad platform interior" gate for the checked seeds.
  - The legacy generated-world province diversity status remains
    `generated_world_gate_incomplete` on all three seeds because it still
    includes passive-margin lowland object acceptance.  That is intentionally
    separated from the P104A internal-mosaic gate.
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py aevum/modules/terrain.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104a_continental_mosaic_expression_fixture -q`
    -> `1 passed in 1.22s`.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104a_generated_world_continental_mosaic_gate_pass -q`
    -> `1 passed in 114.57s`.

Next:

- Run `P104B.crust_sediment_interior_elevation_repair`: use crust thickness,
  sediment thickness, erosion state, basement age, and province type to reduce
  broad high-flat interiors and to create more Earth-like basins, shields,
  platforms, old uplands, and lowlands.
- Define and execute `P104C.internal_geographic_block_initialization`: large
  continents should contain inherited craton/platform/mobile-belt/basin/terrane
  blocks from initialization and lifecycle evolution, so diversity is not only a
  final-stage rendering repair.
- Keep P105 planform/ribbon/coastline repairs separate from P104B/P104C unless
  generated-world evidence shows a direct coupling.

2026-06-27 - P104B crust/sediment interior elevation repair completed

- Production work:
  - Added a post-P104A P104B surface repair:
    `TerrainModule._apply_crust_sediment_interior_elevation_repair(...)`.
  - The repair consumes production province codes, final continental detail,
    crust thickness, sediment thickness, crust stability, and crust domain.
  - It caps unsupported high-flat platforms, shields, and basins while leaving
    explicit orogen and LIP plateau highlands available for the landform
    inventory.
  - It includes a P103-compatible continental land-floor recovery guard for
    Earth-like worlds.  The guard restores low-lying submerged continental
    cells near the sea-level floor instead of re-raising high interiors.
  - The P104A platform fallback split now also handles platform-dominated
    exposed components with only one-cell interior width.  Fallback relief is
    low-amplitude and deterministic, based on crust/sediment ranking rather
    than random texture.
- Diagnostic work:
  - Added `P104B.high_flat_land_floor_gate`.
  - Added `P104B.mosaic_preservation_gate`.
  - Added `test_p104b_generated_world_interior_elevation_gate_pass`.
  - Updated the P104A fixture to recognize intentional overlarge-platform
    fallback splitting.
- Current output:
  - `PYTHONPATH=. .venv/bin/python -m aevum.diagnostics.tectonics_bench --suite P104B --out out_bench_p104b_crust_sediment_interior_elevation_20260627`
    -> `status: pass`, deterministic `True`.
  - Summary and CSV:
    `out_bench_p104b_crust_sediment_interior_elevation_20260627/tectonics_bench_summary.json`
    and
    `out_bench_p104b_crust_sediment_interior_elevation_20260627/p104b_microbenchmarks.csv`.
- Acceptance evidence:
  - Seeds checked: `42`, `31415`, `16180`; resolution `900` cells.
  - Minimum continental land fraction: `0.2500309`.
  - Maximum direct high-flat interior fraction of continental land: `0.0081642`
    against the `0.020` gate.
  - Minimum basin/lowland fraction of continental land: `0.5112364`.
  - Minimum lowland `<500m` fraction: `0.4890401`.
  - Minimum lowland `<1000m` fraction: `0.6090758`.
  - Maximum unparented highland fraction of highlands: `0.0714688`.
  - Maximum P104B adjusted area fraction: `0.0533709`.
  - Maximum mean absolute P104B surface delta: `125.24 m`.
  - Maximum land-floor recovery area fraction: `0.0110784`.
  - P104A mosaic after P104B remains valid:
    min province classes per major component `4`,
    max largest internal province fraction `0.607173`,
    failing major component count `0`.
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/terrain.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104a_continental_mosaic_expression_fixture -q`
    -> `1 passed in 1.19s`.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104a_generated_world_continental_mosaic_gate_pass -q`
    -> `1 passed in 117.17s`.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104b_generated_world_interior_elevation_gate_pass -q`
    -> `1 passed in 117.00s`.

Next:

- Define `P104C.internal_geographic_block_initialization` in detail and then
  implement it.  The goal is to reduce reliance on late fallback splitting by
  giving large continents inherited internal craton, platform, mobile-belt,
  basin, rift, and terrane blocks from initialization and lifecycle processes.
- Keep checking P104A/P104B together: P104C should improve natural causality
  without reintroducing high-flat interiors, all-platform major components, or
  land-floor regressions.
- Leave planform/ribbon/coastline repairs to P105 unless P104C evidence shows
  direct coupling.

2026-06-27 - P104C internal geographic block initialization completed

- Implementation:
  - Added an inherited continental block layer in `TectonicsModule`:
    `tectonics.internal_geographic_block_id`,
    `tectonics.internal_geographic_block_code`, and
    `tectonics.internal_geographic_blocks`.
  - Block classes now include craton core, stable platform, intracratonic
    basin, mobile belt, rifted margin, and accreted terrane.
  - Initialization is deterministic and physically parameterized from
    proto-crust potential, component width, proto-craton masks, margin
    potential, crustal thickness, and stability.
  - Lifecycle rules inherit blocks with moving crust parcels, clear them on
    new oceanic crust, and rewrite them at collision, extension, and arc/terrane
    accretion zones.
  - Added a maturation rule for overwide old mobile/accreted collages so large
    interiors can settle into covered platforms and intracratonic basins rather
    than remaining one undifferentiated active belt.
  - Terrain now consumes inherited blocks in the production province graph, and
    P104A has a detail-recovery step for subdued old mobile belts, basins,
    rifts, and accreted terranes that would otherwise render as broad platform.
- Evidence:
  - `P104C.internal_geographic_block_gate` passed on seeds `42`, `31415`, and
    `16180` at `900` cells.
  - Minimum internal block object count: `35`.
  - Minimum internal block class count: `5`.
  - Minimum major-continent internal block classes: `4`.
  - Maximum largest internal block class fraction: `0.461525`.
  - Minimum terrain internal-block consumption area fraction: `0.124434`.
  - Minimum detail recovery from internal blocks: `0.0566894`.
  - Maximum late P104A platform fallback split area fraction: `0.0`.
  - P104A mosaic after P104C remained valid: minimum province classes per major
    component `6`, maximum largest internal province fraction `0.379324`, and
    failing major component count `0`.
  - P104B after P104C remained valid: minimum continental land fraction
    `0.250010`, maximum high-flat interior fraction `0.00823075`, minimum
    basin/lowland fraction `0.559750`, and minimum lowland `<500m` fraction
    `0.461995`.
- Outputs:
  - `out_bench_p104c_internal_geographic_block_initialization_20260627/tectonics_bench_summary.json`
  - `out_bench_p104c_internal_geographic_block_initialization_20260627/p104c_microbenchmarks.csv`
  - `out_bench_p104b_after_p104c_regression_20260627/tectonics_bench_summary.json`
  - `out_bench_p104b_after_p104c_regression_20260627/p104b_microbenchmarks.csv`
- Verification:
  - `PYTHONPATH=. .venv/bin/python -m py_compile aevum/modules/tectonics.py aevum/modules/terrain.py aevum/diagnostics/tectonics_bench.py tests/test_tectonics_bench.py aevum/features.py`
    -> pass.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104c_internal_geographic_block_initialization_gate_pass -q`
    -> `1 passed in 121.18s`.
  - `PYTHONPATH=. .venv/bin/pytest tests/test_tectonics_bench.py::test_p104b_generated_world_interior_elevation_gate_pass -q`
    -> `1 passed in 122.42s`.
- Interpretation:
  - The immediate "大陆内部太单一" failure is now causally addressed in
    tectonics rather than only corrected in terrain rendering.
  - The next required evidence is visual: generate higher-resolution worlds and
    compare `elevation.png`, continental detail/province maps, and `hexmap.png`
    for natural-looking internal geographic regions without striping or
    checkerboard artifacts.

2026-06-27 - P104C 2500-cell visual audit recorded

- Assets:
  - `out_p104c_visual_audit_3_worlds_20260627/p104c_visual_contact_sheet.png`
  - `out_p104c_visual_audit_3_worlds_20260627/p104c_visual_audit_summary.json`
  - Per-world folders for seeds `42`, `31415`, and `16180` with
    `elevation.png`, `continental_detail_provinces.png`,
    `internal_geographic_blocks.png`, and `hexmap.png`.
- Metrics:
  - Minimum internal block object count: `76`.
  - Minimum major-continent internal block classes: `6`.
  - Maximum largest internal block class fraction: `0.404946`.
  - Minimum terrain internal-block consumption area fraction: `0.134815`.
  - Minimum detail recovery from internal blocks: `0.0544012`.
  - P104A fallback platform split area fraction: `0.0` for all three worlds.
  - Minimum visible province classes per major component: `5`.
  - Maximum largest internal province fraction: `0.402752`.
  - Compiler envelope passed for all three worlds.
- Review:
  - Continental interiors now have visible craton/platform/basin/mobile-belt/
    rift/accreted-terrain differentiation at `2500` cells.
  - The visual maps still show diagonal striping and cell-scale patchiness in
    internal block/detail expression.  This appears to come from repeated
    parcel rasterization of process cargo, not from a lack of block diversity.
  - A direct truth-cargo smoothing experiment was rejected because it harmed
    seed `31415` high-flat/lowland behavior.  The fix should be a province
    aggregation/render-expression layer or a more physical object grouping,
    not pre-terrain smoothing of `tectonics.internal_geographic_block_code`.
- Remaining evidence gaps:
  - Seed `31415` has lowland `<500m` fraction `0.275115` at `2500` cells.
  - Seed `16180` has continental land fraction `0.238825` at `2500` cells.
  - These do not invalidate P104C's internal-diversity claim, but they prevent
    declaring the full plate/terrain stage closed without further scaling
    checks.

2026-06-27 - P104D internal block regional expression and floor/lowland audit recorded

- Assets:
  - `out_bench_p104d_internal_block_region_expression_20260627/`
  - `out_p104d_visual_audit_3_worlds_floor_lowland_20260627/p104d_visual_contact_sheet.png`
  - `out_p104d_visual_audit_3_worlds_floor_lowland_20260627/p104d_visual_audit_summary.json`
- Implementation:
  - Added `terrain.internal_geographic_block_region_code` as a render-scale
    regional expression of inherited internal blocks.
  - The region field does not mutate `tectonics.internal_geographic_block_code`
    and does not drive elevation, avoiding the rejected P104C smoothing
    regression.
  - Added rendering for both raw internal block cargo and regional expression.
  - Added mature earthlike continental floor restoration for the 2500-cell
    seed `16180` land-floor failure.
  - Added P104B earthlike lowland recovery for the seed `31415` lowland and
    basin/lowland failure.
- 900-cell benchmark evidence:
  - `P104D` passed `2/2` microbenchmarks, deterministic `True`.
  - Minimum same-neighbor improvement: `0.209270`.
  - Maximum regional tiny-block area fraction: `0.00383257`.
  - Minimum regional class count: `5`.
  - P104B/P104C preservation remained valid.
- 2500-cell visual-audit evidence across seeds `42`, `31415`, and `16180`:
  - Minimum same-neighbor improvement: `0.209081`.
  - Maximum regional tiny-block area fraction: `0.00285908`.
  - Minimum regional block class count: `6`.
  - Maximum largest regional block fraction: `0.487588`.
  - Minimum continental land fraction: `0.250021`.
  - Maximum high-flat interior fraction: `0.00738038`.
  - Minimum basin/lowland fraction: `0.415050`.
  - Minimum lowland `<500m` fraction: `0.364727`.
  - Compiler envelope passed for all three worlds.
- Interpretation:
  - Continental interiors now have enough inherited block diversity at object
    and regional-expression scale for the checked generated worlds.
  - The remaining visual defect is in `terrain.continental_detail`: detail
    provinces still show diagonal cell-cargo striping, even when the regional
    internal block layer is coherent.
- Next:
  - P104E/P105 should make detail-province and compiler-visible terrain
    expression consume region/province objects rather than raw per-cell cargo.

2026-06-27 - P104E continental detail regional expression audit recorded

- Assets:
  - `out_bench_p104e_continental_detail_region_expression_20260627/`
  - `out_p104e_visual_audit_3_worlds_detail_region_20260627/p104e_visual_contact_sheet.png`
  - `out_p104e_visual_audit_3_worlds_detail_region_20260627/p104e_visual_audit_summary.json`
- Implementation:
  - Added `terrain.continental_detail_region_code` as a compiler-visible
    regional expression of continental detail provinces.
  - The field is derived from raw `terrain.continental_detail`, regional
    internal geographic blocks, and production continental province codes.
  - Raw detail remains exported separately as
    `continental_detail_raw_provinces.png`; the standard
    `continental_detail_provinces.png` now shows the regional expression.
  - The strategy hex compiler now consumes the regional detail field, avoiding
    direct propagation of raw per-cell cargo striping into `hexmap.png`.
- Formal 900-cell evidence:
  - `P104E` passed `2/2` microbenchmarks, deterministic `True`.
  - Minimum same-detail-neighbor improvement: `0.261236`.
  - Maximum regional tiny-detail area fraction: `0.0355363`, down from raw
    maximum `0.268735`.
  - Minimum detail-region class count: `7`.
  - Maximum detail-region largest class fraction: `0.482486`.
  - Compiler and P104B/P104D preservation checks passed.
- 2500-cell visual-audit evidence across seeds `42`, `31415`, and `16180`:
  - Minimum same-detail-neighbor improvement: `0.287927`.
  - Maximum regional tiny-detail area fraction: `0.0236155`, down from raw
    maximum `0.189791`.
  - Minimum detail-region class count: `7`.
  - Maximum largest detail-region class fraction: `0.691492`.
  - Minimum continental land fraction: `0.250021`.
  - Maximum high-flat interior fraction: `0.00738038`.
  - Minimum basin/lowland fraction: `0.415050`.
  - Minimum lowland `<500m` fraction: `0.364727`.
  - Compiler envelope passed for all three worlds.
- Interpretation against real-Earth geomorphology goals:
  - The display/compiler layer now has coherent continental provinces instead
    of diagonal cell-scale artifacts, so the previous "striped detail map"
    failure is mostly addressed.
  - This does not yet fully solve continental interior geomorphology.  Real
    continents have multiple physiographic provinces within platforms and
    basins: shields, sedimentary basins, old eroded orogens, escarpments,
    rift shoulders, plateaus, forelands, and inherited uplands.  P104E
    classifies those provinces more coherently, but elevation still lacks a
    sufficiently rich province-scale response in some large interiors.
- Next research target:
  - Define and benchmark a province-scale elevation/landform response step that
    uses regional detail, terrain provinces, basin-fill age, inherited uplands,
    rift shoulders, foreland loads, passive-margin escarpments, and denudation
    state.  The goal is to reduce large single-color interior surfaces while
    preserving realistic lowland/platform/basin dominance.

2026-06-27 - P104F inland geomorphology elevation-region response audit recorded

- Assets:
  - `out_bench_p104f_after_old_orogen_detail_split_20260627/`
  - `out_bench_p104e_after_p104f_old_orogen_split_20260627/`
  - `out_p104f_visual_audit_3_worlds_old_orogen_split_20260627/p104f_visual_contact_sheet.png`
  - `out_p104f_visual_audit_3_worlds_old_orogen_split_20260627/p104f_visual_audit_summary.json`
- Implementation:
  - Added `terrain.inland_geomorphology_region_code` as a production field
    that expresses inland physiographic regions on exposed continental crust.
  - Added a bounded province-scale surface response after P104B/P104E-style
    regionalization.  The response consumes regional continental detail,
    continental province code, crust thickness/stability, sediment, rift
    potential, inherited uplands, and current elevation.
  - Recomputed `terrain.continental_detail_region_code` after the inland
    geomorphology response, allowing over-broad old-orogen detail regions to
    split into platform, basin, rift-basin, shield, and orogen expression where
    the inland region state supports that split.
  - Added rendering for `inland_geomorphology_regions.png` and a P104F
    generated-world microbenchmark suite.
- Formal 900-cell evidence:
  - `P104F` passed `2/2` microbenchmarks, deterministic `True`.
  - Minimum inland geomorphology region class count: `3`.
  - Minimum applied continental inland fraction: `0.896415`.
  - Maximum adjusted surface area fraction: `0.0731802`.
  - Maximum mean absolute surface delta: `40.4078 m`.
  - Minimum major-continent elevation-band count after P104F: `3`.
  - Maximum largest major-continent elevation-band fraction after P104F:
    `0.487019`.
  - Minimum major-continent relief IQR after P104F: `361.746 m`.
  - P104B/P104E preservation checks passed, including compiler consistency,
    land-mask preservation, lowland/platform floors, and high-flat suppression.
- P104E regression after old-orogen split:
  - `P104E` still passed `2/2` microbenchmarks, deterministic `True`.
  - The compiler still consumes `terrain.continental_detail_region_code`.
  - Maximum 900-cell detail-region largest class fraction after the split:
    `0.526635`.
- 2500-cell generated-world visual audit across seeds `42`, `31415`, and
  `16180`:
  - Minimum inland geomorphology region class count: `5`.
  - Minimum applied continental inland fraction: `0.871044`.
  - Maximum adjusted surface area fraction: `0.0787871`.
  - Maximum mean absolute surface delta: `89.9002 m`.
  - Minimum major-continent elevation-band count after P104F: `4`.
  - Maximum largest major-continent elevation-band fraction after P104F:
    `0.388827`.
  - Minimum major-continent relief IQR after P104F: `604.762 m`.
  - Maximum P104E detail-region largest class fraction fell to `0.533063`,
    down from the previous P104E 2500-cell maximum of `0.691492`.
  - Minimum continental land fraction remained `0.250020`, maximum high-flat
    interior fraction remained `0.00736949`, and compiler envelope passed for
    all three worlds.
- Interpretation against real-Earth geomorphology goals:
  - The immediate continental-interior monotony failure is now addressed by a
    process-linked inland region response rather than by simple recoloring.
    Large interiors can now contain lowlands, platforms, basins, old eroded
    uplands, rift-related interiors, and shield-like areas with distinct
    elevation behavior.
  - This is closer to real continents, where most continental area is not
    high mountain terrain but still contains multiple physiographic provinces.
  - The current output should not be promoted as final Earth-realistic terrain
    yet.  Remaining visual residuals include some broad smooth platform/lowland
    surfaces, local diagonal/parcel texture, and limited explicit drainage/
    basin-controlled relief expression.
- Next:
  - Keep climate, ocean circulation, and monsoon redesign paused.
  - If another terrain pass is needed, target drainage/basin and surface-process
    expression inside the existing P104F inland regions rather than adding more
    random detail classes.
  - Re-run a higher-resolution promotion audit before declaring the
    plate/terrain stage closed.
