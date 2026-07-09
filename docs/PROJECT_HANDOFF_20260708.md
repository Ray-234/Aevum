# Aevum Project Handoff - 2026-07-08

This note is the short entry point for continuing Aevum on a CPU cluster.
The long-running design logs remain in the other `docs/*.md` files; this file
only records the current state, what is intentionally frozen, and what should
run next.

## Current State

- The source tree, tests, feature registry, small reference fixtures, and
  development plans are ready to move to GitHub.
- The plate/terrain generation stage is considered good enough to close for
  now.  It is not perfect, but the current priority has moved away from map
  shaping toward climate/biome derivation from accepted terminal terrain.
- Existing generated maps, benchmark outputs, videos, and external raw climate
  datasets are local reproducible artifacts, not repository source.  They are
  ignored to keep the GitHub repository usable.
- The in-repo fast Aevum climate engine is frozen as a production climate path.
  It remains useful as a diagnostic/prototype layer, but it should not be
  further tuned to hide Earth-replay residuals.
- Climate work should proceed by fitting one real-Earth subgraph at a time,
  visually comparing Earth reference maps against replay maps before any metric
  promotion.

## Included In The Repository

- `aevum/`: engine, modules, compiler, renderer, diagnostics.
- `tests/`: regression, release-gate, climate-reference, and replay tests.
- `docs/`: current plans and status archives.
- `data/registry/features.yaml`: feature contracts.
- `data/reference/earth_hypsometry_fixture_20260627.json`: small Earth
  hypsometry fixture.
- `data/reference/earth_climate/source_manifest.json`: manifest of Earth
  climate reference sources gathered for calibration.
- `data/reference/etopo5/`: small topography reference files used during
  Earth-geomorphology comparison.

## Excluded Local Artifacts

The local workspace contains roughly 15 GB of generated and downloaded files.
They are excluded by `.gitignore`, especially:

- `/out*/`: benchmark outputs, rendered PNGs, videos, arrays, replay packets.
- `/data/reference/earth_climate/raw/`: WorldClim, NOAA, OSCAR, GODAS,
  Koppen-Geiger, ecoregion, and related downloaded source files.
- `/data/reference/earth_climate/cache/` and
  `/data/reference/earth_climate/processed/`: derived climate-reference cache.
- `*.npz`, `*.png`, `*.mp4`, and similar generated binary artifacts outside
  explicitly allowed documentation assets.

Do not commit these to normal GitHub history.  Several raw files exceed
GitHub's 100 MB object limit.  If the cluster needs persistent large artifacts,
use a separate object store, Git LFS, or a documented data-sync step.

## Plate And Terrain Closure

The active plate/terrain closure checklist is archived in
`docs/PLATE_TECTONICS_ENGINEERING_PLAN.md`.  The most important accepted
properties are:

- terminal maps are plausibly Earth-like enough for downstream climate tests;
- plate-boundary semantics have object layers for ridges, trenches, transforms,
  convergent parents, orogenic hierarchy, shelves, slopes, abyssal plains,
  island arcs, microcontinents, and ocean-floor objects;
- ridge/trench/bathymetry/crust-age logic is much more process-linked than the
  early noise-driven maps;
- 8000-cell maps are the practical calibration target, with 24000-cell maps
  used for visual promotion and 72000-cell selected snapshots for high-detail
  experiments.

Known residual terrain issues are deferred unless they block climate:

- mountain spines and branch ranges can still be smoother and more continuous;
- some high-latitude and edge classifications remain overactive;
- grid texture is still visible in some object layers;
- deep-time plate reorganization is still more heuristic than a full physical
  plate model.

The current instruction for climate experiments is: use the accepted terminal
terrain as input and do not keep reshaping plate generation unless a climate
failure is clearly caused by terrain.

## Climate Direction

The current climate conclusion is that temperature, precipitation, winds, and
ocean currents cannot be fitted independently by a quick rule engine.  The
correct route is:

1. Use real Earth as the calibration target before generated worlds.
2. Fit one subgraph at a time, in dependency order:
   terrain/land-sea geometry -> energy/SST boundary state -> pressure centers
   -> winds -> ocean currents -> SST/sea ice -> moisture transport ->
   precipitation -> Koppen -> biome.
3. For each subgraph, render the real Earth reference map, the Aevum/external
   replay map, and a residual map on the same projection and scale.
4. Read the maps first and assign visible residuals to upstream owners before
   changing code or parameters.
5. Treat global means and aggregate scores as regression guards, not as the
   primary proof of correctness.

Detailed plans:

- `docs/EARTH_BASED_CLIMATE_FITTING_PLAN.md`
- `docs/CLIMATE_MECHANISM_MODELING_PLAN.md`
- `docs/CLIMATE_COUPLING_RESEARCH_NOTES.md`
- `docs/EARTH_CLIMATE_REFERENCE_CALIBRATION_PLAN.md`

## External Climate Engine Gate

A separate local sandbox under `/Users/rayw/claude/geoworld` has been used to
test whether a real GCM-like tool should replace the frozen fast climate
engine.  That sandbox is not part of this repository, but its current lesson is
important for Aevum:

- ExoPlaSim is usable as the first external climate-engine gate.
- The zero-obliquity copied setup was rejected because it removed seasons.
- The corrected Earth-orbit setup uses Earth-like obliquity, eccentricity,
  solar constant, CO2, mixed-layer depth, and sea ice.
- A 10-year T21 Earth cold-start run completed, but it was still drifting cold
  and was not a pass.  The next work is a default-physics audit and controlled
  parameter/preset sweep around stable PlaSim defaults.
- MPAS remains a possible later engine, but it is heavier to configure.  Do not
  jump to MPAS before the ExoPlaSim/Koppen/biome postprocessing contract is
  proven end to end.

For generated Aevum worlds, the intended production path is eventually:

```text
accepted terminal terrain + land/sea + orbital/atmospheric parameters
  -> external climate engine climate normals
  -> Koppen postprocessing from 12 monthly T/P fields
  -> BIOME4-style or DGVM-lite biome inference
  -> optional vegetation/climate iteration
```

## Local Smoke Commands

Install and run tests:

```bash
python3 -m venv .venv
./.venv/bin/pip install -e .
./.venv/bin/python -m pytest -q
```

Local verification note from this handoff:

- On this macOS workspace, plain Python 3.13 test collection initially failed
  because `pyexpat` loaded the system `libexpat`.  The working invocation was:

```bash
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib \
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/opt/expat/lib \
./.venv/bin/python -m pytest -q
```

- Follow-up fix, 2026-07-09: the three observed failures above were repaired.
  The repair also fixed a later `tests/test_engine.py` seam-audit false
  positive where polar longitude-wrap edges were counted as antimeridian
  dateline seams.
- Verification after the fix:
  `tests/test_engine.py::test_parented_hotspot_chain_can_form_limited_archipelago`,
  `tests/test_engine.py::test_tidally_locked_world_does_not_get_earthlike_seasons`,
  and
  `tests/test_engine.py::test_ocean_currents_are_basin_constrained_and_transport_heat`
  passed together (`3 passed in 54.04s`).
- Related regression subset passed:
  `tests/test_core.py`, `tests/test_p173_ocean_lifecycle_gate.py`,
  `tests/test_earth_climate_ocean_spatial_gate.py`, and
  `tests/test_earth_climate_coupled_consistency_gate.py`
  (`22 passed in 1.03s`).
- Seam regression checks passed:
  `tests/test_engine.py::test_tectonic_diagnostics_cover_current_world_and_archive`
  and `tests/test_engine.py::test_truth_layers_are_cyclic_at_dateline`
  (`2 passed in 123.63s`).
- A full test-suite pass has still not been completed locally because the
  suite is long-running on this workstation.  The cluster should rerun the full
  suite before large parameter sweeps.

Generate a quick world:

```bash
./.venv/bin/python -m aevum.cli run \
  --preset earthlike \
  --cells 8000 \
  --frames 36 \
  --out out_smoke_earthlike_8000
```

Run the six-world release-style gate:

```bash
./.venv/bin/python -m aevum.cli p12 \
  --presets earthlike waterworld arid stagnant_lid tidally_locked frozen \
  --cells 8000 \
  --frames 36 \
  --out out_p12_cluster_smoke
```

Profile resolutions before large sweeps:

```bash
./.venv/bin/python -m aevum.cli profile-resolution \
  --preset earthlike \
  --cells 8000 24000 \
  --frames 36 \
  --out out_profile_cluster
```

## Immediate Next Step On The Cluster

Start with reproducibility and throughput:

1. Clone the repository and create a clean virtual environment.
2. Run the full test suite once.
3. Run `profile-resolution` at 8000 and 24000 cells.
4. Run a small process-level parallel seed sweep for accepted terrain presets.
5. Only after the terrain baseline is reproducible, start external climate
   engine experiments and Earth-subgraph replay fitting.

Use process-level parallelism.  Do not share writable output directories across
workers.
