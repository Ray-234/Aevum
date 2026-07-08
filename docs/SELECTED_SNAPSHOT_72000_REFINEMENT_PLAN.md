# Selected Snapshot 72000 Refinement Plan

Status: active.

This plan covers high-resolution refinement for one selected terminal snapshot.
It is intentionally separate from deep-time plate generation: 8000-24000 cell
runs own plate topology, major continent/ocean layout, crustal parentage,
orogen belts, trenches, ridges, basins, and lifecycle diagnostics.  The 72000+
layer must inherit those parent fields and add local physiographic detail
without rewriting the parent topology.

## Goal

Generate a 72000-cell refinement package from a chosen P107 terminal array
archive.  The package must make selected maps visually richer at local scale
while preserving the parent snapshot's major land/ocean, basin, plate, crust,
orogen, trench, ridge, and semantic object structure.

## Hard Rules

- Do not run a full 4500 Myr history at 72000 cells for routine iteration.
- Use a 24000-cell terminal array archive as the parent source unless a lower
  resolution parent is explicitly selected for a smoke test.
- Preserve parent land/ocean topology by default.  Any sign flips must be
  process-limited and reported separately.
- Process-backed micro-island or atoll-islet promotion is opt-in only.  The
  default refinement remains topology-preserving; opt-in promotion must write
  rank, mask, delta, and sign-flip metrics so it can be accepted or rejected
  independently.
- Store parent elevation, refined elevation, and detail delta together.
- Every visual iteration must include image inspection of elevation,
  bathymetry, parent/refined comparison, and detail delta.
- Texture phases may be deterministic, but amplitudes must be controlled by
  inherited process fields: orogen hierarchy, continental detail, crust age,
  ocean depth province, boundary masks, and object masks.

## Initial Pipeline

1. Load a P107 terminal metrics JSON and its `p107_terminal_arrays.npz`.
2. Reconstruct the parent `SphereGrid` and create a 72000-cell target grid.
3. Resample parent fields to target cells:
   - continuous fields use inverse-distance interpolation;
   - discrete fields and masks use nearest-parent assignment;
   - elevation interpolation preserves the nearest parent land/ocean sign.
4. Derive a process-conditioned `detail_delta_m`:
   - orogen hierarchy/spine/halo/apron adds narrow high-relief texture;
   - continental shield/platform/basin/rift/orogen/plateau codes set inland
     relief amplitudes;
   - young ocean crust, ridge masks, trench masks, seamount chains, oceanic
     plateaus, microcontinents, and island arcs set bathymetric detail
     amplitudes.
5. Clamp refined elevation to preserve parent sign, except for explicitly
   reported process-backed island candidates.
6. Write a P107-compatible refined array archive and metrics JSON.
7. Render the normal P107 array QA pack plus 72000-specific comparison maps.
8. Record diagnostics: land fraction change, sign-flip fraction, detail
   percentiles, parent source path, and process mask area fractions.

## Visual Review Checklist

- Continents keep the same first-order shape as the parent.
- Mountain belts remain aligned with parent orogen hierarchy and do not become
  isolated speckle.
- Lowlands and basins gain subtle relief without becoming high plateaus.
- Shelves, slopes, abyssal plains, ridges, trenches, plateaus, and seamount
  chains remain distinguishable in bathymetry.
- Added details reduce blockiness rather than adding diagonal grid artifacts.
- Delta maps show detail concentrated in process-backed regions, not uniform
  global noise.

## Progress Log

- 2026-07-06: Started selected-snapshot refinement after full-history 8000
  preflight showed roughly 10 minute runtime for one terminal frame and is not
  a practical 72000 iteration path.  Chosen route: derive 72000 cells from an
  existing 24000 P107 terminal array parent, then inspect generated images
  before further parameter work.
- 2026-07-06: Added first terminal microgeomorphology layer for selected
  snapshots: land coast-distance, flow accumulation, river rank, conservative
  visible river path rank, small lake masks, and river-mouth delta masks.  The
  accepted default output is
  `out_selected_snapshot_72000_refinement_seed707_v13_default_river_objects_20260706`.
  It preserves land/ocean sign, runs in about 20 seconds for 72000 cells, and
  removes the earlier dotted/comb-like river rendering by drawing sparse
  main-stem and tributary path objects.
- 2026-07-06: Rejected an experimental basin-outlet A* river pass.  It produced
  longer lowland rivers, but visual review showed looped/U-shaped edge paths and
  runtime increased from about 20 seconds to 64-126 seconds per 72000 run.  Do
  not re-enable this as a default patch.  The long-river problem should be
  solved later with explicit drainage-basin/outlet objects, not ad hoc global
  path search.
- 2026-07-06: Added conservative marine microgeomorphology.  The accepted output
  is `out_selected_snapshot_72000_refinement_seed707_v18_marine_class_colorbar_20260706`.
  It writes `field__selected_snapshot_marine_delta_m`,
  `field__selected_snapshot_ocean_coast_distance_passes`,
  `field__selected_snapshot_shelf_break_rank`,
  `mask__selected_snapshot_reef_atoll`, and
  the marine shoal masks split into union, seamount, oceanic plateau,
  microcontinent, and island-arc classes.  It also adds three QA renders:
  `selected_snapshot_marine_microgeomorphology.png`,
  `selected_snapshot_marine_zoom_sheet.png`, and
  `selected_snapshot_marine_object_classes.png`.
- 2026-07-06: Rejected the first marine attempt
  `out_selected_snapshot_72000_refinement_seed707_v14_marine_micro_20260706`.
  It overcorrected broad ocean depth: ocean detail p95 rose above 1000 m,
  shelf-break overlays covered roughly one fifth of ocean cells, and mid-ocean
  ridge fallback was mislabeled as seamount shoals.  The accepted v18 pass keeps
  land/ocean sign flips at zero, lowers ocean detail p95 to about 498 m, keeps
  marine nonzero p95 near 348 m, and restricts shoals to actual seamount,
  oceanic plateau, microcontinent, and island-arc objects outside the polar
  edge zone.
- 2026-07-06: Added explicit drainage-basin objects without reintroducing A*
  routing.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v24_smooth_trunk_paths_20260706`.
  It writes `field__selected_snapshot_drainage_basin_id`,
  `field__selected_snapshot_basin_trunk_rank`, and
  `field__selected_snapshot_floodplain_rank`, plus
  `selected_snapshot_drainage_basins.png`.  The first basin attempt,
  `out_selected_snapshot_72000_refinement_seed707_v19_drainage_basins_20260706`,
  improved river lengths but produced long polar/edge trunk lines.  The accepted
  v24 pass restricts basin objects outside the polar edge zone and increases
  spacing between selected trunks.  Compared with v18, visible river path cells
  rise from 323 to 494, path maximum coast-distance rises from 8 to 16 passes,
  and trunk p95 coast-distance is about 13 passes while land/ocean sign flips
  remain zero.
- 2026-07-06: Added a conservative meander/floodplain presentation layer.  It
  writes `field__selected_snapshot_meander_belt_rank` and uses that rank to
  smooth selected basin-trunk paths during rendering.  This is intentionally a
  display/microgeomorphology pass: endpoints remain on the receiver tree and no
  new basin routing is introduced.  The accepted v24 hydrology zoom sheet is
  smoother than v21 while retaining the same terrain, coastline, and basin
  metrics.
- 2026-07-06: Added explicit lake-basin and river-mouth morphology.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v25_lake_delta_morphology_20260706`.
  It writes `field__selected_snapshot_lake_basin_rank`,
  `field__selected_snapshot_delta_fan_rank`, and
  `mask__selected_snapshot_delta_plain`.  Visual review shows sparse lakes now
  have local basin halos and river mouths have process-bounded delta fan/plain
  expression without changing parent land/ocean topology.  Metrics remain
  conservative: land/ocean sign flips stay at zero, lake-basin cells cover about
  0.46% of land, delta fans about 0.78% of ocean, and delta plains about 1.72%
  of land.
- 2026-07-06: Added process-backed island and atoll candidate QA without
  changing parent land/ocean topology.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v27_island_atoll_candidate_latgate_20260706`.
  It writes `field__selected_snapshot_island_candidate_rank`,
  `field__selected_snapshot_atoll_candidate_rank`,
  `mask__selected_snapshot_process_island_candidate`, and
  `mask__selected_snapshot_atoll_candidate`, plus
  `selected_snapshot_island_atoll_candidates.png`.  The first candidate pass,
  `out_selected_snapshot_72000_refinement_seed707_v26_island_atoll_candidates_20260706`,
  was rejected because its strongest islet candidate was selected from a
  high-latitude edge seamount belt near 68 degrees.  The accepted v27 pass adds
  a default latitude suitability gate for islet candidates: candidate cores are
  now confined to about -27 to +33 degrees in this seed, while atoll candidates
  remain tropical/subtropical.  Sign flips remain zero.
- 2026-07-06: Added conservative reef/atoll morphology ranks while preserving
  parent topology.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v30_reef_atoll_zoom_20260706`.
  It writes `field__selected_snapshot_reef_rim_rank`,
  `field__selected_snapshot_atoll_lagoon_rank`, and
  `field__selected_snapshot_fringing_reef_rank`, plus
  `selected_snapshot_reef_atoll_morphology.png` and
  `selected_snapshot_reef_atoll_zoom_sheet.png`.  The first morphology pass,
  `out_selected_snapshot_72000_refinement_seed707_v28_reef_atoll_morphology_20260706`,
  was superseded because lagoon cells outnumbered reef-rim cells, which did not
  read like an atoll.  The accepted v30 pass uses a higher-scoring atoll core as
  lagoon and the remaining candidate/halo cells as reef rim.  In this seed,
  reef rim covers about 0.10% of ocean cells, lagoon about 0.03%, and fringing
  reef about 0.07%; land/ocean sign flips remain zero and ocean detail p95 is
  effectively unchanged.
- 2026-07-06: Added topology-preserving coastal morphology ranks.  The accepted
  output is
  `out_selected_snapshot_72000_refinement_seed707_v33_coastal_plain_retune_20260706`.
  It writes `field__selected_snapshot_coastal_delta_m`,
  `field__selected_snapshot_coastal_plain_rank`,
  `field__selected_snapshot_coastal_cliff_rank`,
  `field__selected_snapshot_shoreface_rank`,
  `field__selected_snapshot_barrier_lagoon_rank`, and
  `field__selected_snapshot_estuary_rank`, plus
  `selected_snapshot_coastal_morphology.png` and
  `selected_snapshot_coastal_zoom_sheet.png`.  The first coastal pass,
  `out_selected_snapshot_72000_refinement_seed707_v31_coastal_morphology_20260706`,
  was rejected because coastal plain/cliff and shoreface ranks were too broad
  and leaked into polar/edge coasts.  The accepted v33 pass gates coastal ranks
  by latitude, margin setting, active-margin exclusion, and local relief.  In
  this seed, coastal cliffs cover about 1.6% of land, shoreface about 10.9% of
  ocean, barrier/lagoon about 1.1% of ocean, estuary about 0.8% of ocean, and
  coastal plains remain sparse because the selected world has very little
  low-relief non-orogenic coast.  Land/ocean sign flips remain zero.
- 2026-07-06: Added topology-preserving submarine highland ranks.  The accepted
  output is
  `out_selected_snapshot_72000_refinement_seed707_v35_submarine_highlands_tuned_20260706`.
  It writes `field__selected_snapshot_submarine_highland_delta_m`,
  `field__selected_snapshot_seamount_peak_rank`,
  `field__selected_snapshot_seamount_apron_rank`,
  `field__selected_snapshot_oceanic_plateau_edge_rank`, and
  `field__selected_snapshot_abyssal_hill_field_rank`, plus
  `selected_snapshot_submarine_highlands.png` and
  `selected_snapshot_submarine_highlands_zoom_sheet.png`.  The first pass,
  `out_selected_snapshot_72000_refinement_seed707_v34_submarine_highlands_20260706`,
  was rejected for a small high-latitude microcontinent/plateau-edge leak and a
  too-dense abyssal-hill field.  The accepted v35 pass keeps seamount peaks at
  about 0.20% of ocean cells, seamount aprons about 0.79%, plateau edges about
  1.71%, and abyssal-hill fields about 2.05%; land/ocean sign flips remain zero.
- 2026-07-06: Added an explicit opt-in process-island/atoll-islet promotion
  gate.  The default route still preserves topology; `--allow-process-islands`
  now promotes only high-rank process-backed islet and atoll-rim candidates
  after marine/reef morphology has been computed.  It writes
  `field__selected_snapshot_process_island_promotion_rank`,
  `field__selected_snapshot_process_island_promotion_delta_m`,
  `mask__selected_snapshot_process_island_promoted`, and
  `mask__selected_snapshot_atoll_islet_promoted`, plus
  `selected_snapshot_process_island_promotion.png` and
  `selected_snapshot_process_island_promotion_zoom_sheet.png`.  The first
  opt-in pass
  `out_selected_snapshot_72000_refinement_seed707_v36_process_island_promotion_20260706`
  was superseded because the unified top-N promotion selected no atoll islets.
  The accepted opt-in audit output is
  `out_selected_snapshot_72000_refinement_seed707_v37_process_island_promotion_tuned_20260706`:
  it promotes 15 parent-ocean cells, including 5 atoll-islet cells, with
  land/ocean sign flips at about 0.0208% of all cells and promoted relative
  elevations around +90 to +236 m.  This remains an opt-in selected-snapshot
  promotion gate, not the default topology-preserving baseline.
- 2026-07-06: Added a subcell island/atoll microshape QA layer.  The accepted
  output is
  `out_selected_snapshot_72000_refinement_seed707_v38_island_atoll_microshapes_20260706`.
  It writes `field__selected_snapshot_islet_microshape_rank` and
  `field__selected_snapshot_atoll_microshape_rank`, plus
  `selected_snapshot_island_atoll_microshapes.png` and
  `selected_snapshot_island_atoll_microshapes_zoom_sheet.png`.  This pass does
  not change elevation or topology relative to the v37 opt-in promotion audit:
  sign flips remain about 0.0208%, promoted parent-ocean cells remain 15, and
  atoll-islet promoted cells remain 5.  The new QA layer draws small elongated
  island-chain symbols and reef-ring/lagoon symbols in zoomed views so the
  selected-snapshot output no longer needs to visually read as only coarse
  full-cell island blocks.
- 2026-07-06: Added coastal process linework for the river-mouth/shelf
  transition.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v40_coastal_linework_tuned_20260706`.
  It writes `field__selected_snapshot_delta_distributary_rank`,
  `field__selected_snapshot_estuary_funnel_rank`, and
  `field__selected_snapshot_barrier_spit_rank`, plus
  `selected_snapshot_coastal_process_linework.png` and
  `selected_snapshot_coastal_process_linework_zoom_sheet.png`.  The first
  linework pass
  `out_selected_snapshot_72000_refinement_seed707_v39_coastal_process_linework_20260706`
  was superseded because zoomed river-mouth symbols drew too many adjacent
  rank cells and read as arrow clusters.  The accepted v40 pass keeps the same
  terrain and rank coverage, but renders only locally spaced process peaks:
  delta distributary and estuary-funnel ranks cover about 0.77% of ocean cells,
  barrier-spit ranks about 0.87%, and land/ocean sign flips remain unchanged at
  about 0.0208% in the opt-in micro-island audit.
- 2026-07-06: Added fluvial/lacustrine microshape QA for inland rivers and
  lakes.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v42_fluvial_lacustrine_microshapes_tuned_20260706`.
  It writes `field__selected_snapshot_meander_scroll_rank`,
  `field__selected_snapshot_floodplain_swale_rank`, and
  `field__selected_snapshot_lake_shoreline_rank`, plus
  `selected_snapshot_fluvial_lacustrine_microshapes.png` and
  `selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png`.  The first
  pass
  `out_selected_snapshot_72000_refinement_seed707_v41_fluvial_lacustrine_microshapes_20260706`
  was superseded because meander-scroll ranks were too low for the render
  threshold and did not visibly draw.  The accepted v42 pass keeps terrain and
  topology unchanged while making low-amplitude scroll bars, floodplain swales,
  and lake shorelines visible in zoomed QA.  In this selected seed, scroll
  ranks cover about 0.23% of land, swales about 2.60%, lake shorelines about
  0.25%, and opt-in sign flips remain about 0.0208%.
- 2026-07-06: Started carrying the fluvial/lacustrine QA semantics into the
  refined elevation array as low-amplitude selected-snapshot microrelief.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v45_fluvial_microrelief_tuned_20260706`.
  It writes `field__selected_snapshot_fluvial_microrelief_delta_m` and renders
  `selected_snapshot_fluvial_microrelief_delta.png` plus
  `selected_snapshot_fluvial_microrelief_zoom_sheet.png`.  The first v43 pass
  was superseded because the microrelief attribution compared the final clamped
  terrain against an unclamped baseline, making sea-level clamps look like
  tens-of-meters microrelief.  The v44 attribution fix was visually safe but
  too conservative.  The accepted v45 pass keeps the same parent topology and
  opt-in sign flips as v42, applies actual lowland river/lake microrelief to
  about 2.91% of land cells, and keeps the applied nonzero p95 at about 7 m
  with a local range of roughly -19 m to +7 m.
- 2026-07-06: Started carrying coastal process linework into the refined
  elevation array as topology-preserving shallow-water / river-mouth
  microrelief.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v46_coastal_process_microrelief_20260706`.
  It writes `field__selected_snapshot_coastal_process_microrelief_delta_m` and
  renders `selected_snapshot_coastal_process_microrelief_delta.png` plus
  `selected_snapshot_coastal_process_microrelief_zoom_sheet.png`.  This pass
  keeps the same parent topology and opt-in sign flips as v45, applies coastal
  process microrelief to about 1.58% of ocean cells, and keeps the applied
  nonzero p95 at about 32 m with a local range of roughly -15 m to +35 m.
  Delta/estuary cells primarily show shallow river-mouth channel incision;
  barrier-spit/lagoon cells show local shoal buildup.
- 2026-07-07: Retuned the river-mouth coastal process split after visual QA
  showed that delta-distributary and estuary-funnel ranks were aliases of the
  same `delta_fan_rank` cells.  The current verification output is
  `out_selected_snapshot_72000_refinement_seed707_v69_coastal_estuary_delta_split_20260707`.
  Delta support now favours open, shallow shelf mouths; estuary support favours
  more confined or embayed coast cells.  The v68 intermediate fixed the severe
  sign error where delta mouths were cut downward; v69 reduces delta/estuary
  overlap from effectively all high-score cells to 67 cells, with 177
  delta-only and 115 estuary-only high-score cells in this seed.  Coverage
  stays conservative: delta-distributary cells cover about 0.34% of ocean,
  estuary funnels about 0.25%, barrier spits about 0.61%, coastal process
  microrelief about 1.04%, and nonzero process abs-p95 is about 21 m.  Visual
  review shows separate delta and estuary zoom windows, positive delta-mouth
  deposition, and local estuary incision.  Remaining weakness: these are still
  sparse cell-scale symbols and deltas do not yet form strong trunk-to-mouth
  distributary networks.
- 2026-07-07: Re-rendered the v69 coastal QA sheets with hydrology-aware
  coastal process linework.  The renderer now uses saved `river_receiver`,
  `river_path_rank`, basin-trunk rank, and flow accumulation to draw only
  short river-mouth trunk hints whose downstream land cell touches a
  delta/estuary coastal process object.  It remains presentation-only and does
  not change the 72000 arrays or parent topology.  The first visual attempt was
  too dominant and painted long blue river strokes through barrier examples;
  the accepted retune halves the near-mouth trace length, reduces opacity and
  linewidth, and suppresses mouths where barrier-spit support dominates
  river-mouth support.  Visual review now shows delta/estuary panels with a
  readable trunk-to-mouth cue while keeping barrier-spit panels focused on
  coastal shoal objects.  Remaining weakness: this is still a QA linework hint,
  not a generated subcell delta-plain network with avulsion lobes, tidal
  channels, and distributary-mouth bars.
- 2026-07-07: Added actual river-mouth subfeatures to the coastal process
  arrays.  The current verification output is
  `out_selected_snapshot_72000_refinement_seed707_v71_coastal_mouth_bar_tidal_channel_tuned_20260707`.
  It writes `field__selected_snapshot_delta_mouth_bar_rank` and
  `field__selected_snapshot_estuary_tidal_channel_rank`, then folds those
  ranks into `field__selected_snapshot_coastal_process_microrelief_delta_m`.
  Mouth-bar rank is derived from open-shelf delta support and adds weak
  positive shoal relief; tidal-channel rank is derived from confined estuary
  support and adds weak negative channel relief.  The v70 first attempt was
  visually too strong because mouth bars stacked on top of the existing delta
  process delta and frequently hit the +30 m cap.  The accepted v71 retune
  lowers the added amplitudes; in seed 707 the new ranks cover about 0.15% and
  0.12% of current-ocean cells respectively, the process nonzero abs-p95 is
  about 25 m, and land/ocean sign flips stay at about 0.0208%.  Visual review
  shows mouth-bar and tidal-channel semantics in both the microrelief raster
  and coastal process linework without broadening the coastal footprint.
  Remaining weakness: this is still cell-scale microrelief, not a true
  subcell delta-plain network with avulsion lobes, channel bifurcation,
  tidal flats, and mouth-bar complexes.
- 2026-07-06: Started carrying island/atoll subcell QA symbols into the
  refined elevation array as sparse opt-in island-chain / reef-shoal
  microrelief.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v48_island_atoll_microrelief_tuned_20260706`.
  It writes `field__selected_snapshot_island_atoll_microrelief_delta_m` and
  renders `selected_snapshot_island_atoll_microrelief_delta.png` plus
  `selected_snapshot_island_atoll_microrelief_zoom_sheet.png`.  The first v47
  pass was superseded because almost every nonzero cell hit the +54 m cap and
  read as hard isolated uplift points.  The accepted v48 pass lowers the cap,
  adds a one-cell weak halo, keeps the same parent topology and opt-in sign
  flips as v46, applies island/atoll microrelief to about 0.23% of ocean cells,
  and keeps the applied nonzero p95 at about 38 m.  In this seed the added
  contribution is mostly positive shoal/islet buildup; lagoon deepening is
  already represented by the earlier reef/atoll morphology delta.
- 2026-07-06: Added a sparse shelf-edge / upper-slope microrelief pass after
  opt-in process-island promotion and before island/atoll microrelief.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v50_shelf_slope_microrelief_tuned_20260706`.
  It writes `field__selected_snapshot_shelf_slope_microrelief_delta_m` and
  renders `selected_snapshot_shelf_slope_microrelief_delta.png` plus
  `selected_snapshot_shelf_slope_microrelief_zoom_sheet.png`.  The first v49
  pass was rejected because it painted a too-broad, mostly positive shelf-edge
  band across about 8.66% of ocean cells and risked recreating shallow marine
  patches.  The accepted v50 pass tightens the shelf-break/depth gate, raises
  the seed percentile, weakens the halo, and keeps applied relief small: about
  4.24% of ocean cells receive a nonzero contribution, nonzero p95 is about
  3.3 m, and the applied range is roughly -6.1 m to +6.0 m after topology
  clamps.  Visual review shows thin shelf-edge and upper-slope line segments
  rather than new broad shoal objects.
- 2026-07-06: Added a first conservative deep-ocean fabric pass for selected
  transform/fracture and abyssal-plain settings.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v53_deep_ocean_fabric_line_tuned_20260706`.
  It writes `field__selected_snapshot_deep_ocean_fabric_delta_m`,
  `field__selected_snapshot_fracture_zone_rank`, and
  `field__selected_snapshot_abyssal_plain_fabric_rank`, plus
  `selected_snapshot_deep_ocean_fabric_delta.png` and
  `selected_snapshot_deep_ocean_fabric_zoom_sheet.png`.  The v51 pass was
  rejected because transform/fracture masks rendered as broad blob-like
  patches; v52 reduced that too far into sparse bead-like ranks.  The accepted
  v53 pass keeps the real elevation contribution very small and topology-safe:
  about 5.16% of ocean cells receive a nonzero contribution, nonzero p95 is
  about 4.9 m, fracture-zone rank covers about 0.83% of ocean cells, abyssal
  plain fabric rank covers about 4.34%, and land/ocean sign flips remain
  unchanged.  Remaining weakness: this is a safe raster fabric layer, not a
  true continuous fracture-zone polyline system.  The latter should be repaired
  upstream in boundary/object linework.
- 2026-07-07: Retuned the deep-ocean fracture fabric after v71 visual review
  showed that fracture-zone ranks still appeared as short beads and red/blue
  rings in zoomed QA.  The current verification output is
  `out_selected_snapshot_72000_refinement_seed707_v76_deep_ocean_fracture_trough_only_20260707`.
  The accepted code extracts line-like fracture candidates, fills only
  single-cell gaps, removes the positive fracture shoulder uplift, and adds a
  QA-only principal-axis line overlay for fracture components.  The v73/v74
  attempts were rejected because they bridged too aggressively and turned
  fracture fabric into broad blobs.  The accepted v76 path is a middle ground:
  fracture-rank connected components drop from 158 to 91, tiny components
  (`<=3` cells) drop from 125 to 34, fracture-rank coverage rises moderately
  from about 0.83% to about 1.25% of ocean cells, deep-ocean fabric nonzero
  abs-p95 stays low at about 5.5 m, and land/ocean sign flips remain unchanged
  at about 0.0208%.  Visual review shows more readable fracture-trough axes
  without making the actual bathymetry a broad shoal/trough patch.  Remaining
  weakness: the source fracture objects are still raster masks, so this is a
  gap-filled 72000 fabric cue, not a true boundary-derived fracture-zone
  polyline network.
- 2026-07-07: Retuned hydrology and fluvial/lacustrine QA rendering after v76
  visual review showed that the `river mouth` zoom could center on an offshore
  delta-fan peak with the river pushed to the crop edge, and that lake
  shorelines still rendered as scattered or mesh-like symbols.  The current
  verification output reuses
  `out_selected_snapshot_72000_refinement_seed707_v76_deep_ocean_fracture_trough_only_20260707`
  as a render-only QA package.  The accepted code now scores river-mouth zoom
  centers from landward coastal outlet cells with adjacent offshore delta-fan
  support, applies a temperate-latitude preference to zoom selection, reduces
  lake-basin river overlay clutter, and draws lake-shoreline linework only
  along cells adjacent to the actual lake mask.  This is a presentation/QA
  improvement only: saved terrain, lake, river, and delta arrays are unchanged.
  Regression tests now assert that river-mouth zooms choose the landward mouth
  and lake-shoreline zooms prefer a representative lower-latitude candidate
  over a stronger high-latitude one.
- 2026-07-07: Retuned island/atoll candidate semantics after v76/v77 visual
  review and distance telemetry showed that atoll, reef-rim, and micro-island
  candidates were still dominated by nearshore / embayed shallow seas.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v79_offshore_atoll_islet_rebalance_20260707`.
  The generator now separates offshore island/atoll candidates from nearshore
  fringing reefs: island candidates require process-backed offshore support,
  atoll candidates require a stricter open-ocean gate, and hybrid rim/island
  promotion can still classify a sparse opt-in atoll islet when reef-rim
  support is clear.  In the selected seed, atoll candidates drop from about
  0.108% to about 0.042% of ocean cells, but their nearshore (`d<=2`) count
  drops from 43 cells to zero; atoll lagoons drop from 16 cells to 6 cells and
  all remaining lagoon cells are `d>=5`; fringing reefs remain nearshore as
  intended.  Process-island promotion drops from 16 cells to 8 cells, all
  `d>=5`, with one sparse atoll islet retained.  Land/ocean sign flips improve
  from about 0.0208% to about 0.0083%.  Visual review shows the island/atoll
  QA sheets now read as offshore volcanic-highland / reef candidates rather
  than inner-bay stickers.  Remaining weakness: ring-lagoon geometry is still
  symbolic and cell-scale; true reef-rim arcs and island coastlines still
  belong to a later local-super-resolution pass.
- 2026-07-07: Retuned shelf-edge / upper-slope microrelief after v79/v76 visual
  review showed that the layer still read as dotted terrace/channel speckles
  rather than continuous shelf-break and upper-slope process axes.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v80_shelf_axis_continuity_20260707`.
  The generator now derives a shelf-slope axis rank from shelf-break support,
  bridges supported one-cell gaps, and derives low-amplitude terrace/channel
  deltas from that axis instead of directly from isolated texture peaks.  The
  ocean-cell footprint rises only moderately from about 4.28% to about 5.26%,
  while nonzero p95 amplitude drops from about 3.32 m to about 2.26 m and
  land/ocean sign flips remain unchanged at about 0.0083%.  Component telemetry
  confirms less dotting: at a 0.75 m threshold, small components (`<=3` cells)
  drop from 235 to 185, while the largest connected component grows from 14 to
  31 cells.  Visual review shows more readable shelf-break line segments and
  upper-slope channels without painting continuous broad shoals along every
  coastline.  Remaining weakness: this is still a 72000-cell axis/rank proxy,
  not explicit canyon polylines or sediment-routing bathymetry.
- 2026-07-07: Retuned submarine seamount-chain expression after v80 visual
  review showed that seamount peaks still rendered as small filled orange
  patches rather than discrete volcano/guyot-like highs arranged along a chain.
  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v81_seamount_chain_axis_20260707`.
  The generator now spaces high-rank seamount peak seeds, bridges those seeds
  into a weak chain axis, and derives the apron from the chain axis while
  keeping peak rank core-only.  Seamount peak cells drop from 98 to 37 and are
  all single-cell peak cores; seamount apron cells remain broad enough to carry
  the chain context, dropping only from 397 to 308 with a max component still
  at 48 cells.  Submarine-highland nonzero p95 amplitude drops from about
  20.24 m to about 17.28 m, and land/ocean sign flips remain unchanged at
  about 0.0083%.  Visual review shows seamount examples as bead-like peaks
  embedded in a low-amplitude apron chain rather than one saturated blob.
  Remaining weakness: guyot flat tops, explicit volcanic-line polylines, and
  detailed plateau escarpment planforms still require local super-resolution or
  upstream object-line geometry.
- 2026-07-07: Retuned oceanic plateau / microcontinent escarpment expression
  after v81 visual review showed that plateau edges could still read as
  complete geometric frames.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v82_plateau_escarpment_segments_20260707`.
  The generator now selects texture-supported escarpment segments from the
  plateau/microcontinent boundary instead of forcing the full object outline
  into the rank field.  Plateau-edge footprint drops from about 1.68% to about
  1.04% of ocean cells; strong edge cells (`rank>0.75`) drop from 345 to 54,
  and the largest strong component drops from 65 to 16 cells.  Submarine-
  highland nonzero p95 amplitude drops from about 17.28 m to about 16.02 m,
  and land/ocean sign flips remain unchanged at about 0.0083%.  Visual review
  shows the former closed-frame edge replaced by shorter escarpment segments
  while retaining weak plateau/microcontinent context.  Remaining weakness:
  the segments are still cell-scale rank proxies, not explicit plateau
  escarpment polylines or guyot/flat-top morphology.
- 2026-07-06: Added a first lowland / alluvial microrelief pass on the land
  side, after fluvial/lacustrine microrelief and before coastal morphology.
  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v55_lowland_alluvial_microrelief_tuned_20260706`.
  It writes `field__selected_snapshot_lowland_alluvial_microrelief_delta_m`,
  `field__selected_snapshot_alluvial_fan_rank`,
  `field__selected_snapshot_lowland_plain_rank`, and
  `field__selected_snapshot_piedmont_apron_rank`, plus
  `selected_snapshot_lowland_alluvial_microrelief_delta.png` and
  `selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png`.  The v54
  pass was rejected as too broad because platform-like lowlands could enter the
  plain rank without enough basin/floodplain/lake/delta support, covering about
  28.4% of land.  The accepted v55 pass gates the layer to mid/low latitudes
  and process-supported lowlands: nonzero lowland/alluvial microrelief covers
  about 16.0% of land cells, nonzero p95 is about 9.4 m, alluvial fan ranks
  cover about 0.80% of land, lowland plain ranks about 9.3%, and piedmont
  apron ranks about 7.0%.  Land/ocean sign flips remain unchanged.  Remaining
  weakness: alluvial fans are still rank patches, not physically routed fan
  planforms; that should be handled by later source-to-sink sediment routing or
  local super-resolution rather than increasing this delta.
- 2026-07-06: Added a first land-side coastal depositional microrelief pass
  after shoreface/barrier/estuary morphology.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v57_coastal_depositional_microrelief_narrow_20260706`.
  It writes
  `field__selected_snapshot_coastal_depositional_microrelief_delta_m`,
  `field__selected_snapshot_coastal_depositional_plain_rank`,
  `field__selected_snapshot_strandplain_rank`, and
  `field__selected_snapshot_tidal_flat_rank`, plus
  `selected_snapshot_coastal_depositional_microrelief_delta.png` and
  `selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png`.  The
  first v56 pass was rejected as semantically too broad: it was low-amplitude,
  but ordinary shoreface-adjacent lowlands could trigger the object and about
  15.4% of land cells received a nonzero contribution.  The accepted v57 pass
  removes generic shoreface-only triggering and requires lowland support,
  river-mouth/delta support, or barrier/lagoon support.  Coverage drops to
  about 4.1% of land cells, nonzero p95 is about 5.0 m, depositional plain
  ranks cover about 4.2% of land, strandplain ranks about 2.6%, and tidal-flat
  ranks about 1.9%.  Land/ocean sign flips remain unchanged.  Remaining
  weakness: these are still cell-scale depositional hints; subcell beach-ridge,
  tidal-channel, and marsh/lagoon geometry belongs in the later local
  super-resolution pass.
- 2026-07-06: Retuned submarine highland rank semantics after visual review of
  the v57 zoom sheet.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v58_submarine_highland_rank_tuned_20260706`.
  No new arrays were added; this pass tightens
  `field__selected_snapshot_oceanic_plateau_edge_rank` and
  `field__selected_snapshot_abyssal_hill_field_rank`.  The specific bug fixed
  here was that microcontinent cells could be forced into the plateau-edge
  rank across their interior rather than only along object margins, causing
  blocky filled patches.  The v58 pass only gives high edge rank to actual
  microcontinent/plateau boundary cells, while allowing a weak one-cell halo.
  Abyssal-hill fields now use a stricter seed threshold plus a small halo, so
  they read as local clusters instead of isolated single-cell speckles.  Metrics
  remain conservative: plateau-edge rank is about 1.68% of ocean cells,
  abyssal-hill rank about 2.36%, submarine-highland nonzero p95 remains about
  20.2 m, and land/ocean sign flips are unchanged.  A regression test now
  asserts that high microcontinent edge rank does not fill object interiors.
- 2026-07-06: Retuned lake-shoreline microrelief after visual review showed
  lake shoreline symbols and the applied fluvial/lacustrine delta could still
  read as isolated lake-core dots rather than a local shoreline ring.  The
  accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v59_lake_shoreline_ring_tuned_20260706`.
  The change keeps the same fields and assets, but changes
  `field__selected_snapshot_lake_shoreline_rank` so it is seeded from the
  outer lake rim, excludes `mask__selected_snapshot_lakes`, and applies only a
  low-amplitude positive shore-step delta.  Metrics stay conservative:
  lake-shoreline rank increases from about 0.25% to about 0.37% of land cells,
  fluvial/lacustrine nonzero p95 falls slightly to about 6.3 m, and land/ocean
  sign flips are unchanged.  A regression test now asserts that lake shoreline
  rank does not cover lake-core cells.  Remaining weakness: lakes are still
  terminal-snapshot local masks, not true subcell waterbody polygons with
  routed inflow/outflow shore geometry.
- 2026-07-06: Retuned hydrology QA linework rendering after visual review of
  v58/v60 showed that weak river segments still appeared as short blue dashes
  on the full hydrology map and in lake-basin zooms.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v61_hydrology_linework_render_tuned_20260706`.
  No refined array fields were added; this is a presentation/QA tune for
  `selected_snapshot_hydrology.png`,
  `selected_snapshot_hydrology_zoom_sheet.png`, and fluvial/lacustrine overlay
  sheets.  The renderer now draws `river_path_rank` as connected paths with
  minimum path length, uses stricter thresholds for full-map weak river
  segments, and applies local spacing to avoid dense short-line clutter.  The
  v59/v60 hydrology-linework outputs were intermediate checks; v61 keeps zoomed
  major rivers, lake basins, and river mouths readable while making the full
  hydrology map less dash-like.  Metrics and topology remain effectively
  unchanged: land/ocean sign flips remain about 0.000208 and the hydrology
  fields are not regenerated by this renderer-only pass.  Remaining weakness:
  lake-basin zooms can still show multiple true basin-trunk paths crossing a
  view.  Fixing that requires basin/outlet hierarchy or local super-resolution,
  not looser render filtering.
- 2026-07-06: Retuned reef/atoll QA overlays after visual review showed that
  reef rim, atoll lagoon, and fringing-reef ranks were hard to distinguish
  against shallow-water elevation colors.  The accepted output is
  `out_selected_snapshot_72000_refinement_seed707_v63_reef_atoll_overlay_points_20260706`.
  No refined array fields changed; v62 was an intermediate fixed-color overlay
  pass, and v63 adds sparse rank point anchors over the same fixed semantic
  colors.  `selected_snapshot_reef_atoll_morphology.png` and
  `selected_snapshot_reef_atoll_zoom_sheet.png` now use gold circles for reef
  rims, blue circles for lagoon cells, and green triangles for fringing reefs,
  while keeping the underlying rank rasters translucent.  Reef/atoll metrics
  are unchanged from v61: reef rim about 0.10% of ocean cells, lagoon about
  0.03%, fringing reef about 0.07%, and land/ocean sign flips remain about
  0.000208 with opt-in process-island promotion enabled.  Remaining weakness:
  this is still QA symbology over 72000-cell ranks, not subcell reef-rim or
  lagoon coastline geometry.
- 2026-07-06: Retuned island/atoll microrelief amplitude after visual review of
  v63/v64 showed that promoted islet and atoll cells still rendered as saturated
  red blobs in `selected_snapshot_island_atoll_microrelief_zoom_sheet.png`.
  The current focused output is
  `out_selected_snapshot_72000_refinement_seed707_v65_island_atoll_microrelief_small_core_20260706`.
  The generator now keeps halo response weak, gates reef/fringing buildup by
  the island/atoll object axis, and caps this late microrelief contribution at
  -18..20 m.  Process-island promotion remains responsible for any opt-in
  emergence; island/atoll microrelief is now only a small selected-snapshot
  texture layer.  Metrics preserve topology and object counts: land/ocean sign
  flips remain about 0.000208, reef rim about 0.10% of ocean cells, lagoon
  about 0.03%, and fringing reef about 0.07%.  The island/atoll microrelief
  nonzero p50 drops from about 3.02 m in v63 to about 0.47 m in v65, and p95/max
  drop from 38 m to 20 m.  The full CLI render path still has a slow tail; v65
  arrays and core metrics were produced by the CLI, while the deep-ocean,
  submarine-highland, island/atoll, and reef/atoll QA images were rendered back
  from the saved v65 arrays.
- 2026-07-07: Added grouped selected-snapshot QA rendering and a render-only
  recovery path for already generated refinement arrays.  `selected-snapshot-refine`
  now accepts `--render-groups`, with groups such as `base`, `hydrology`,
  `marine`, `shelf`, `deep-ocean`, `submarine`, `island-atoll`, and `coastal`;
  non-`all` runs can also skip the p107 contact-sheet render.  The new
  `selected-snapshot-render-groups` CLI reads an existing
  `selected_snapshot_refinement_metrics.json` / `selected_snapshot_refined_arrays.npz`
  pair and renders requested QA groups without recomputing the 72000 refinement.
  The verification output
  `out_selected_snapshot_72000_refinement_seed707_v66_grouped_island_submarine_20260706`
  was interrupted after arrays/metrics were written, then recovered through the
  render-only CLI.  It now has 11 `submarine` + `island-atoll` refinement
  assets, no p107 assets, and a `selected_snapshot_render_group_summary.json`.
  The submarine-highland zoom-center heuristic was also retuned so the seamount
  panel prefers offshore, still-submerged peaks instead of the highest nearshore
  shoal rank; visual review now shows the seamount example in a deep-ocean
  setting.

## Current Visual Assessment

- 72000 refinement now preserves parent continent/ocean topology and improves
  local land/ocean texture without global stripe artifacts.
- Hydrology is useful as a first selected-snapshot micro layer: lakes and deltas
  are sparse, rivers are object-like instead of dense raster noise, and
  hydrology deltas do not dominate the terrain field.  The v61 QA linework
  retune also reduces full-map short-dash clutter by drawing connected path
  objects and filtering weak isolated segments.  The v77 render-only QA retune
  improves the inspection sheets: river-mouth crops now center on landward
  outlets rather than offshore fan peaks, and lake shorelines render as local
  rim linework rather than scattered dots or triangular mesh.
- Marine microgeomorphology is now useful as a first conservative pass: shelves
  and shelf breaks read more clearly, reefs/atolls are sparse tropical shallow
  objects, and seamount/plateau shoals no longer paint continuous mid-ocean
  ridge belts.  Marine object class maps now split seamount, oceanic plateau,
  microcontinent, island-arc, shelf-break, and reef/atoll classes.  The v66
  zoom-center tune makes the submarine-highland QA sheet sample a more realistic
  offshore seamount peak rather than a nearshore shelf high.
- Hydrology now has a first explicit drainage-basin object layer.  It improves
  long lowland trunk rivers without global path search and avoids the rejected
  polar trunk artifacts from v19.  Remaining weakness: trunks are still
  receiver-tree paths and local tributary geometry can still look angular in
  zoomed views, but selected lowland trunks now render as smoothed meander-belt
  paths.  A later pass should add true basin outlet hierarchy and floodplain
  corridor widths rather than changing parent topology.
- Lake and delta objects are now explicit enough for selected-snapshot QA:
  lake-basin halos are local rather than global wetland noise, and delta fans
  remain tied to high-accumulation coastal outlets.  Remaining weakness: small
  islands and atolls are still primarily represented as shallow marine objects.
  The next safe step is to expose process-backed island/atoll candidates as
  explicit QA objects while preserving parent land/ocean topology by default.
- Process-backed island/atoll candidates are now explicit QA objects.  They are
  deliberately conservative: the default pass ranks candidate islets and atolls
  but does not flip ocean cells into land.  Remaining weakness: candidate
  objects are not yet used to create final selected-snapshot micro-islands,
  lagoons, or fringing-reef coast geometry.  That should remain opt-in until a
  separate island/atoll promotion gate proves that topology changes are local
  and visually useful.
- Reef/atoll morphology now has explicit local ranks for reef rims, atoll
  lagoons, and fringing reefs.  These ranks add small bathymetric relief while
  preserving ocean sign, so they can be audited before any island-promotion
  gate is enabled.  The v63 QA overlay makes the three semantics visually
  distinct enough for review.  The v79 offshore-gate retune separates atoll
  candidates from nearshore fringing reefs: reef rims and lagoons now come from
  offshore process-backed highs, while fringing reefs remain tied to the coast.
  Remaining weakness: the shapes are still cell-scale at 72000 cells and should
  only become visible topographic landforms in a later 72000+ or opt-in
  selected-snapshot promotion pass.
- Coastal morphology now has explicit ranks for coastal plains, coastal cliffs,
  shorefaces, barrier/lagoon shoals, and estuaries.  The accepted pass is
  conservative and topology-preserving: it makes active coasts, shoreface bands,
  and river-mouth estuaries auditable without painting every coastline.  The
  selected seed still lacks broad coastal-plain expression because most low
  coasts are process-classified as active/orogenic or are too narrow at 72000
  cells; this should be revisited with a larger selected snapshot or a
  process-backed depositional/coastal-plain promotion gate.
- Submarine highs now have explicit ranks for seamount peaks, seamount aprons,
  oceanic plateau/microcontinent edges, and abyssal-hill fields.  The v58
  retune is conservative and topology-preserving: plateau/microcontinent edge
  rank no longer fills interiors, and abyssal hills render as sparse local
  clusters instead of isolated single-cell speckles.  The v81 seamount retune
  makes peak cores discrete and shifts continuity into a weak apron/chain axis,
  so seamount examples read more like bead-like volcanic highs along a chain
  instead of filled orange patches.  The v82 plateau retune replaces complete
  outline frames with texture-supported escarpment segments, reducing the
  geometric box impression while keeping weak plateau context.  Remaining
  weakness: guyots, plateau escarpment polylines, and fracture-related volcanic
  ridges are still controlled by 72000-cell object masks rather than true
  subcell bathymetric planforms.
- Shelf-edge and upper-slope microrelief now has a real audited elevation
  contribution rather than only a shelf-break rank map.  The accepted v50 pass
  is intentionally low amplitude and sparse: it adds small terrace/channel
  relief along selected shelf breaks without changing land/ocean topology or
  creating new broad offshore shoals.  The v80 retune makes the same low-energy
  contribution more line-like by deriving a shelf-slope axis rank and bridging
  one-cell gaps, so the zoom sheets read as shelf-edge / upper-slope segments
  rather than isolated dots.  Remaining weakness: this is not a substitute for
  richer deep-ocean object generation or explicit submarine canyon polylines.
  Seamount/guyot chains, fracture zones, abyssal plains, volcanic islands, and
  microcontinental fragments should remain separate semantic/object passes so
  shelf/slope relief does not become a generic ocean texture layer.
- Deep-ocean fabric now has a first explicit selected-snapshot elevation
  contribution and QA layer.  The v53 pass uses parent transform/fault and
  abyssal-plain objects to add very low-amplitude deep-ocean trough/plain
  texture away from ridges, trenches, shelves, and promoted islets.  The v76
  retune reduces fracture-rank fragmentation with one-cell gap filling and
  adds QA linework so fracture-trough axes read as connected features in zoom
  sheets.  Remaining weakness: the source transform/fracture masks are still
  raster blobs or segmented cells; a true continuous fracture-zone polyline
  network still belongs upstream in boundary/object linework rather than in a
  larger 72000 microrelief amplitude.
- Lowland and mountain-front terrain now has a first real selected-snapshot
  microrelief contribution.  The v55 pass derives alluvial fan, lowland plain,
  and piedmont apron ranks from existing basin/floodplain/lake/delta and
  orogenic apron context, then applies a low-amplitude land-only delta.  It
  helps fill the previous gap between river/lake symbols and broad continental
  terrain texture.  Remaining weakness: it does not yet create recognisable
  fan-shaped depositional lobes or fully continuous plains; those require a
  stronger sediment/source-to-sink model or later local super-resolution.
- Land-side coastal depositional terrain now has a first audited microrelief
  contribution.  The v57 pass keeps the accepted shoreface/barrier/estuary
  objects as triggers, but prevents ordinary shoreface adjacency from painting
  all low coasts.  Visual review shows sparse meter-scale delta/strandplain/
  tidal-flat hints tied to low coastal outlets and barrier/lagoon settings.
  Remaining weakness: this does not yet draw continuous beach ridges, tidal
  channels, marsh polygons, or lagoon shoreline geometry; those need local
  super-resolution rather than a broader 72000-cell delta.
- Process-backed micro-island and atoll-islet promotion is now testable without
  changing the default contract.  In the v37 opt-in audit, promoted cells remain
  sparse, low/mid-latitude, and tied to existing shallow marine object ranks.
  The v38 microshape QA layer now expresses those promoted/candidate cells as
  subcell island-chain ellipses and reef-ring/lagoon symbols in rendered zoom
  sheets.  The v48 pass now gives those same symbols a sparse real terrain
  contribution, so island/atoll microshape expression is no longer presentation
  only.  The v65 retune makes that real contribution smaller and less blob-like:
  the object symbols carry subcell shape, while the 72000-cell array carries
  only weak meters-to-tens-of-meters texture.  The v79 retune removes the
  previous nearshore/embayment bias from promoted island and atoll examples:
  process-island promotion is sparser, offshore, and keeps one opt-in atoll
  islet instead of filling shallow coastal bays.  Remaining weakness: this is
  still cell-scale uplift/shoal relief, not
  physical subcell coastline/ring-lagoon geometry, so the promotion gate should
  stay opt-in until a later local super-resolution pass can carry the same
  shapes into a real high-resolution coastal/elevation product.
- River-mouth and coastal-shelf expression now has a first linework QA layer:
  delta distributaries, estuary funnels, and barrier spits are derived from the
  existing hydrology/coastal ranks and rendered as sparse local process peaks.
  The v46 pass now gives those same process ranks a small real bathymetric
  contribution, so river-mouth and barrier/lagoon process expression is no
  longer presentation only.  The v69 retune splits open-shelf delta support
  from confined estuary-funnel support and fixes the previous sign problem that
  made delta mouths read as downcut channels.  The v69 render retune adds
  hydrology-aware near-mouth trunk hints to the coastal process QA sheets.  The
  v71 pass adds actual mouth-bar and tidal-channel ranks to the 72000 arrays.
  Remaining weakness: this is still cell-scale shallow-water microrelief and
  presentation linework rather than explicit subcell coastline / sandbar /
  tidal-channel geometry; a later pass should pair the same semantics with a
  climate/runoff model and local super-resolution coastline geometry.
- Inland hydrology now has a first fluvial/lacustrine microshape QA layer:
  meander scrolls, floodplain swales, and lake shoreline symbols are derived
  from existing basin-trunk, meander-belt, floodplain, and lake-basin ranks.
  The v45 pass now gives those same semantics a small real elevation
  contribution, so selected river/lake microrelief is no longer presentation
  only.  The v59 tune moves lake-shoreline rank off lake cores and onto a
  local outer-rim band, making lake-shore microrelief read more like a shore
  step than isolated lake dots.  The v77 render retune connects only the
  lake-adjacent shoreline rim in QA sheets and adds zoom-center latitude
  preference, so representative lake samples are less likely to be high-polar
  artifacts.  Remaining weakness: this is still cell-scale terminal-snapshot
  microrelief rather than full hydro-geomorphic simulation; receiver-tree
  rivers can still be angular in places, and a later outlet-hierarchy /
  local-super-resolution pass should turn selected river corridors and lake
  basins into explicit subcell planform geometry.
