# CPU Cluster Experiment Plan - 2026-07-08

Target environment: CPU cluster with about 55 available cores.

Primary goal: reproduce the accepted terrain baseline, measure scaling, then
use external climate-engine gates and real-Earth subgraph fitting to derive
Koppen climate maps and biome maps from final Aevum terrain.

## Rules

- Use process-level parallelism only.  Each seed, resolution, or parameter
  preset runs in its own process and output directory.
- Keep deterministic inputs explicit: preset, seed if overridden, cell count,
  frame count, parameter overrides, code commit, and command line.
- Do not use shared mutable caches during a sweep unless the cache is
  read-only for all workers.
- Reserve a few cores for system overhead and I/O.  A practical starting cap is
  48 simultaneous CPU-bound workers on a 55-core node.
- Store generated outputs outside Git history.  Commit only source, plans,
  small fixtures, summaries, and manifest files.

## Stage 0 - Repository And Environment

Acceptance:

- `python -m pytest -q` passes in a clean checkout.
- `aevum.cli presets` lists the six baseline worlds.
- A one-world `run` command writes expected images and JSON outputs.

Suggested commands:

```bash
git clone https://github.com/Ray-234/Aevum.git
cd Aevum
python3 -m venv .venv
./.venv/bin/pip install -U pip
./.venv/bin/pip install -e .
./.venv/bin/python -m pytest -q
./.venv/bin/python -m aevum.cli presets
```

## Stage 1 - Terrain Reproducibility

Run one known-good Earth-like baseline and inspect:

```bash
./.venv/bin/python -m aevum.cli run \
  --preset earthlike \
  --cells 8000 \
  --frames 36 \
  --out out_cluster_stage1_earthlike_8000
```

Acceptance:

- `elevation.png`, `plates.png`, `crust_age.png`, `tectonic_objects.png`,
  `history.png`, and `timeline.png` are generated.
- Terrain is recognizable as the current post-plate-closure quality level:
  coherent landmasses, process-linked mountains, readable ocean floor, no
  obvious seam break.
- Validation checks printed by the command pass.

## Stage 2 - Resolution Scaling

Run a resolution profile before expensive sweeps:

```bash
./.venv/bin/python -m aevum.cli profile-resolution \
  --preset earthlike \
  --cells 8000 24000 \
  --frames 36 \
  --out out_cluster_stage2_resolution_profile
```

Acceptance:

- Record wall time per resolution.
- Record the projected cost for any 72000-cell selected-snapshot work.
- Confirm that 24000-cell renders improve visual continuity without changing
  qualitative world logic.

## Stage 3 - Multi-Preset Baseline Sweep

Run all six presets once at 8000 cells:

```bash
./.venv/bin/python -m aevum.cli p12 \
  --presets earthlike waterworld arid stagnant_lid tidally_locked frozen \
  --cells 8000 \
  --frames 36 \
  --out out_cluster_stage3_p12_8000
```

Acceptance:

- No hard validation failures.
- Warnings are inspected visually, not treated as automatic failures.
- Waterworld and arid worlds should differ by physical regime, not only color
  or sea-level clipping.

## Stage 4 - Process-Level Parallel Seed Sweep

Use independent directories per process.  Example pattern:

```bash
for preset in earthlike waterworld arid; do
  for seed in 42 101 707 909 1001 2024 9001; do
    out="out_cluster_stage4_${preset}_${seed}_8000"
    AEVUM_SEED="$seed" ./.venv/bin/python -m aevum.cli run \
      --preset "$preset" \
      --cells 8000 \
      --frames 36 \
      --out "$out" &
  done
done
wait
```

Current CLI presets carry fixed seeds.  If direct seed override is needed for
large sweeps, add a small explicit CLI seed parameter before running this stage
at scale; do not patch presets manually between runs.

Acceptance:

- Seed-to-seed variation should include different basin and continent layouts.
- Terrain residuals should be classified as either deferred plate-model issues
  or climate-blocking issues.
- Do not reopen plate work unless a repeated defect directly prevents climate
  derivation.

## Stage 5 - External Climate Engine Gate

The in-repo fast climate engine is frozen for production validation.  Use an
external engine gate before generated-world climate production.

Recommended first gate:

1. Reproduce the current ExoPlaSim Earth-orbit T21 baseline.
2. Run a default-physics audit.
3. Run a small parameter/preset sweep around stable defaults.
4. Produce monthly climate normals.
5. Postprocess to Koppen classes from 12 monthly temperature and precipitation
   fields.

Acceptance:

- The Earth run has plausible seasons, zonal gradients, storm-track structure,
  sea ice, and precipitation belts when compared visually with real maps.
- The model is near equilibrium enough for climate normals; one transient year
  is not acceptable for Koppen.
- Each experiment records namelist/config, commit, command, runtime, and
  diagnostic maps.

## Stage 6 - Real-Earth Subgraph Replay

The fitting order remains:

1. terrain and land/sea mask
2. energy/SST boundary state
3. pressure centers
4. winds
5. ocean currents
6. SST and sea ice
7. moisture transport
8. precipitation
9. Koppen classes
10. biome inference

For each subgraph:

- Render the real-Earth reference field.
- Render the model/replay field on the same grid and color scale.
- Render the residual/error field.
- Write a map-read attribution note before code or parameter changes.
- Fix the upstream owner of each visible residual.

Acceptance:

- The visual pattern is right before metrics are used for promotion.
- Metrics can reject regressions but cannot promote a map that is visually
  wrong.
- Downstream fields are observer-only until upstream fields pass.

## Stage 7 - Generated-World Climate And Biomes

Only after the Earth gate is acceptable:

1. Export accepted Aevum terminal terrain and land/sea masks.
2. Run the external climate engine for each generated world.
3. Average enough years for monthly climate normals.
4. Derive Koppen climate classes.
5. Derive biomes through a BIOME4-style model or a documented proxy.
6. Optionally iterate vegetation and climate until the map is stable.

Acceptance:

- Koppen maps are derived from monthly T/P fields, not directly painted from
  latitude or terrain.
- Biome maps are not just Koppen recolors unless explicitly marked as proxy.
- Generated-world climates remain physically consistent with their land/sea,
  topography, orbit, and atmosphere.

## Outputs To Preserve

For every promoted run, keep a small summary bundle:

- command line and environment;
- config/namelist/spec JSON;
- wall time and core count;
- compact diagnostic JSON;
- a contact sheet of key maps;
- map-read attribution note when applicable.

Large raw arrays, model restart files, and videos should live in external
storage or a release artifact, not normal Git history.
