# Aevum Documentation Index

This index is the recommended reading order for new contributors and cluster
experiments.  Most files are research logs; start with the short handoff files
before reading the long plans.

## Start Here

1. [Project handoff](PROJECT_HANDOFF_20260708.md)
   Current state, what is frozen, what is excluded from Git, known test debt,
   and the immediate next steps.
2. [CPU cluster experiment plan](CLUSTER_EXPERIMENT_PLAN_20260708.md)
   How to reproduce terrain baselines, run process-level parallel sweeps, and
   stage external climate-engine work on a 55-core cluster.
3. [Result showcase](RESULT_SHOWCASE.md)
   Curated images and video tracked in Git, with source output directories.
4. [English README](../README_EN.md) and [Chinese README](../README.md)
   Public-facing project summary, quick start, and architecture overview.

## Plate, Terrain, And Geomorphology

- [Plate tectonics engineering plan](PLATE_TECTONICS_ENGINEERING_PLAN.md)
  Main P-series development log and current plate/terrain closure checklist.
- [Plate tectonics refactor plan](PLATE_TECTONICS_REFACTOR_PLAN.md)
  Longer-horizon plan for reducing random heuristics and moving toward
  principle-model tectonics with microbenchmarks.
- [Earth geomorphology coverage](EARTH_GEOMORPHOLOGY_COVERAGE.md)
  Inventory of real-Earth landforms and whether Aevum currently expresses them.
- [Real-Earth geomorphology research plan](REAL_EARTH_GEOMORPHOLOGY_RESEARCH_PLAN.md)
  Research/source/microbenchmark plan for comparing generated worlds to Earth.
- [Continental physiographic architecture plan](CONTINENTAL_PHYSIOGRAPHIC_ARCHITECTURE_PLAN.md)
  Plan for continent-scale physiographic provinces, interiors, basins, plains,
  cratons, mountains, and inherited structure.
- [Historical geomorphology lifecycle plan](HISTORICAL_GEOMORPHOLOGY_LIFECYCLE_PLAN.md)
  How terminal landforms should become process-time objects and archive entries.
- [Selected snapshot 72000 refinement plan](SELECTED_SNAPSHOT_72000_REFINEMENT_PLAN.md)
  High-detail selected-snapshot refinement and visual promotion plan.
- [Tectonics randomness audit](TECTONICS_RANDOMNESS_AUDIT.md)
  Audit of random or overly heuristic parts of the tectonic system.

## Climate And Biomes

- [Earth-based climate fitting plan](EARTH_BASED_CLIMATE_FITTING_PLAN.md)
  Current authoritative workflow: fit one real-Earth subgraph at a time with
  direct visual map comparison and residual attribution.
- [Climate mechanism modeling plan](CLIMATE_MECHANISM_MODELING_PLAN.md)
  Mechanism contracts for terrain, land-sea geometry, pressure, wind, currents,
  SST, moisture, precipitation, Koppen, and biome.
- [Climate coupling research notes](CLIMATE_COUPLING_RESEARCH_NOTES.md)
  Notes on the coupling between topography, wind, ocean currents, temperature,
  and precipitation.
- [Climate system plan](CLIMATE_SYSTEM_PLAN.md)
  Older and broader C-series climate development log.
- [Earth climate reference calibration plan](EARTH_CLIMATE_REFERENCE_CALIBRATION_PLAN.md)
  Earth reference source collection and calibration dataset plan.
- [R2A M0/M1 map-read attribution](R2A_M0_M1_MAP_READ_ATTRIBUTION_20260707.md)
  A concrete example of the map-read attribution workflow for real-Earth
  pressure-source replay.

## Current Execution Guidance

- Treat plate/terrain as the accepted upstream input unless a downstream
  failure is clearly terrain-owned.
- Do not tune the frozen fast climate engine into the final climate product.
- Use real-Earth subgraph replay and external climate engines before generated
  world climate promotion.
- Use process-level parallelism for seed, resolution, and parameter sweeps.
- Keep large generated outputs, raw climate data, and videos outside normal Git
  history unless they are curated showcase assets.
