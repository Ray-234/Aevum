# Tectonics Randomness Audit

Status: R0 audit implemented; P20 plate-split RNG debt removed
Created: 2026-06-22
Latest run: 2026-06-23

This audit records where the tectonics stack still uses randomness.  It is a
working debt ledger for the principle-model tectonics refactor, not a claim that
the current implementation is already physically complete.

Generated audit command:

```text
.venv/bin/python -m aevum.diagnostics.tectonics_bench --suite R0 --out out_bench_r0_p20_lifecycle_final_20260623
```

Generated artifacts:

- `/Users/rayw/Projects/aevum/out_bench_r0_p20_lifecycle_final_20260623/tectonics_bench_summary.json`
- `/Users/rayw/Projects/aevum/out_bench_r0_p20_lifecycle_final_20260623/randomness_inventory.csv`

## Current Counts

R0 currently scans:

- `aevum/modules/tectonics.py`
- `aevum/modules/terrain.py`

Current result:

- total RNG calls: `28`
- allowed: `8`
- temporary: `14`
- forbidden end-state debt: `6`
- unclassified: `0`
- deterministic summary: `true`

The R0 smoke test requires unclassified calls to remain at zero.  Forbidden
calls are allowed only as explicit known debt until the corresponding refactor
stage replaces them.

## Category Meaning

- `allowed`: deterministic RNG stream setup, sub-grid perturbation, bounded
  texture, or representative archive sampling that does not decide physical
  state.
- `temporary`: current implementation debt accepted only until a process model
  replaces it.
- `forbidden`: major geological decision that must not remain random in the
  final principle-model path.
- `unclassified`: missing audit policy; must be classified before further
  tectonics refactor work.

## Forbidden End-State Debt

Initial continent and plate nuclei:

- `_separated_seed_cells`: chooses the first continent nucleus with
  `rng.choice` and perturbs farthest-point selection with `rng.random`.
- `_plate_seed_cells`: chooses plate seeds with random first seed and random
  tie perturbation.
- Replacement: R1/R5 should provide deterministic extrema from mantle heat,
  lithospheric weakness, craton growth, and spacing constraints.

Initial plate motion:

- `TectonicsModule.init_state`: random unit vectors currently assign plate
  Euler poles; random uniform values set initial angular rates.
- Replacement: R2 should compute initial and updated plate motion from torque
  proxies: slab pull, ridge push, basal drag, collision resistance, and mantle
  flow.

## Resolved Former Forbidden Debt

- `TectonicsModule._plume_activity`: formerly selected plume locations with
  `rng.choice`.
- R1 replacement implemented 2026-06-22: plume heads are now selected from
  reproducible local maxima of `tectonics.plume_potential` / mantle upwelling
  fields with deterministic spacing.  R0 no longer reports plume placement as
  a forbidden RNG call.
- `TectonicsModule._refresh_plate_motions`: formerly used
  `random_unit_vectors` and `rng.normal` for pole jitter and rate perturbation.
- R2 replacement implemented 2026-06-23: reorganization refresh now computes
  deterministic torque-proxy pole/rate updates from slab pull, ridge push,
  collision locking, basal drag, and transform friction, with finite motion
  memory.  R0 no longer reports refresh jitter as a forbidden RNG call.
- `TectonicsModule._split_large_plates`: plate split seed, new pole jitter, and
  new rate formerly used RNG.
- P20 replacement implemented 2026-06-23: large-plate split now chooses parent
  plates by resolved-topology pressure, chooses split seeds from rift,
  boundary, crust-type, and stability scores, grows deterministic compact child
  regions, and derives inherited child pole/rate from parent motion and
  parent-child centroid geometry.  R0 no longer reports plate reorganization
  split as a forbidden RNG call.
- P20 capture update implemented 2026-06-23: tiny plate capture now records
  topology/cargo metadata without adding RNG.  Plate labels can be captured
  while continental cargo remains represented by crust and continent/terrane
  ids.
- P20 lifecycle update implemented 2026-06-23: continent and terrane split,
  merge, and capture lineage objects are derived from deterministic id-overlap
  and capture metadata.  R0 counts are unchanged.

## Temporary Debt

Current temporary uses include:

- bounded lognormal continent size targets;
- stochastic continent/craton edge texture;
- legacy stochastic connected-region growth;
- proto-craton extent and thickening variation;
- initial crustal age and stability perturbations that are too important to
  remain arbitrary long-term.

Replacement should occur through R1 thermal fields, R5 continent/margin
lifecycles, and R7 parameter calibration.

## Allowed Uses

Current allowed uses include:

- creating deterministic RNG streams from `RNGKey`;
- small crustal perturbations that can later become sub-grid texture;
- bounded trench recycling age texture;
- representative archive event sampling after physical process masks are
  already determined.

Allowed uses still need guardrails: they must not create, destroy, or place
major geological objects.

## Next Actions

1. Keep `R0.randomness_inventory` in the benchmark suite.
2. Keep `R1` field and plume-trigger microbenchmarks in the benchmark suite.
3. Keep `R2` torque-proxy microbenchmarks in the benchmark suite.
4. Add P20 topology-aware merge/capture and continent lifecycle events so
   topology objects govern both breakup and amalgamation.
5. Then replace continent/plate seed placement with R5-style inherited
   weakness, craton-growth, and deterministic spacing fields.
6. Replace initial Euler pole/rate generation with deterministic mantle-flow,
   slab-pull, ridge-push, and collision-resistance initial conditions.
