"""Replay climate and biome on frozen terminal terrain arrays.

This is the calibration path used after plate/terrain generation is accepted.
It rebuilds a minimal ``WorldState`` from a saved
``terminal_climate_arrays.npz`` archive, restores stellar/orbital forcing from
the matching preset, and reruns only climate plus static biome post-processing.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from aevum import render
from aevum.core.grid import SphereGrid
from aevum.core.module import StepResult
from aevum.core.rng import RNGKey
from aevum.core.state import WorldState
from aevum.core.units import CONSTANTS
from aevum.diagnostics.terminal_climate_biome import (
    _static_npp,
    _summarize_world,
    _write_arrays,
)
from aevum.modules.biosphere import BiosphereModule
from aevum.modules.climate import ClimateModule
from aevum.modules.stellar import StellarModule
from aevum.spec.presets import PRESETS, get_preset


SCHEMA = "aevum.terminal_climate_replay.v1"
FROZEN_INPUT_FIELDS = (
    "terrain.elevation_m",
    "crust.type",
    "tectonics.plate_id",
)


@dataclass(frozen=True)
class TerminalClimateReplayConfig:
    terminal_summary_json: Path
    outdir: Path
    labels: tuple[str, ...] = ()
    render_assets: bool = True


def _json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _apply_step_result(world: WorldState, result: StepResult) -> None:
    delta = result.state_delta
    for key, value in delta.get("fields", {}).items():
        world.fields[key] = np.asarray(value)
    for key, value in delta.get("globals", {}).items():
        world.globals[key] = float(value)
    for key, value in delta.get("networks", {}).items():
        world.networks[key] = value
    for key, value in delta.get("objects", {}).items():
        world.objects[key] = value


def _preset_key_from_summary(summary: dict[str, Any]) -> str:
    raw = str(summary.get("preset", ""))
    if raw in PRESETS:
        return raw
    for key, builder in PRESETS.items():
        if builder().name == raw:
            return key
    label = str(summary.get("assets_dir", "")) + " " + str(summary.get("label", ""))
    for key in PRESETS:
        if key in label:
            return key
    raise ValueError(f"cannot infer preset key from terminal summary preset={raw!r}")


def _world_from_terminal_arrays(summary: dict[str, Any]) -> WorldState:
    arrays_path = Path(summary["arrays"])
    with np.load(arrays_path, allow_pickle=False) as z:
        n = int(np.asarray(z["lat"]).shape[0])
        preset_key = _preset_key_from_summary(summary)
        spec = get_preset(preset_key)
        spec.seed = int(summary.get("seed", spec.seed))
        spec.grid_cells = n
        spec.t_end_myr = float(summary.get("time_myr", spec.t_end_myr))
        grid = SphereGrid.fibonacci(n, spec.radius_m)
        world = WorldState(
            grid=grid,
            spec=spec,
            time_myr=float(summary.get("time_myr", spec.t_end_myr)),
        )
        sea_level = float(np.asarray(z["sea_level_m"], dtype=np.float64).ravel()[0])
        world.set_g("ocean.sea_level_m", sea_level)
        for field in FROZEN_INPUT_FIELDS:
            key = field.replace(".", "__")
            if key in z:
                world.fields[field] = np.asarray(z[key])
        if "terrain.elevation_m" not in world.fields:
            raise ValueError(f"{arrays_path} does not contain terrain.elevation_m")
    return world


def _restore_stellar_forcing(world: WorldState, seed: int) -> dict[str, Any]:
    module = StellarModule()
    key = RNGKey(seed, "terminal_climate_replay_stellar", world.time_myr, 0)
    module.init_state(world, key)
    result = module.step(world, world.time_myr, 0.0, {"terminal_replay": True}, key)
    _apply_step_result(world, result)
    return result.diagnostics


def _run_replayed_climate(world: WorldState, seed: int) -> dict[str, Any]:
    climate = ClimateModule()
    key = RNGKey(seed, "terminal_climate_replay", world.time_myr, 1)
    climate.init_state(world, key)
    result = climate.step(
        world,
        world.time_myr,
        40.0,
        {"terminal_replay": True},
        key,
    )
    _apply_step_result(world, result)

    npp = _static_npp(world)
    biosphere = BiosphereModule()
    biosphere.init_state(world, RNGKey(seed, "terminal_biome_replay", world.time_myr, 2))
    world.fields["biosphere.npp"] = npp
    world.fields["biosphere.biome"] = biosphere._biomes(world, npp)
    return result.diagnostics


def replay_terminal_climate_job(
    summary: dict[str, Any],
    outdir: Path,
    *,
    render_assets: bool = True,
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    seed = int(summary.get("seed", 0))
    world = _world_from_terminal_arrays(summary)
    stellar_diag = _restore_stellar_forcing(world, seed)
    climate_diag = _run_replayed_climate(world, seed)
    if render_assets:
        render.render_world(world, outdir)
    arrays_path = _write_arrays(world, outdir)
    replay_summary = _summarize_world(SimpleNamespace(world=world), climate_diag, outdir)
    replay_summary["source_arrays"] = str(summary["arrays"])
    replay_summary["source_assets_dir"] = str(summary.get("assets_dir", ""))
    replay_summary["stellar_step_diagnostics"] = stellar_diag
    replay_summary["arrays"] = str(arrays_path)
    (outdir / "summary.json").write_text(
        json.dumps(replay_summary, indent=2, default=_json_default)
    )
    return replay_summary


def run_terminal_climate_replay(
    config: TerminalClimateReplayConfig,
) -> dict[str, Any]:
    terminal = json.loads(Path(config.terminal_summary_json).read_text())
    labels = set(config.labels)
    selected: list[dict[str, Any]] = []
    for summary in terminal.get("summaries", []):
        label = Path(str(summary.get("assets_dir", ""))).name
        if labels and label not in labels:
            continue
        item = dict(summary)
        item["label"] = label
        selected.append(item)

    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for summary in selected:
        label = str(summary["label"])
        summaries.append(replay_terminal_climate_job(
            summary,
            outdir / label,
            render_assets=config.render_assets,
        ))
    summaries.sort(key=lambda item: str(item["assets_dir"]))
    batch = {
        "schema": SCHEMA,
        "source_terminal_summary_json": str(config.terminal_summary_json),
        "job_count": len(summaries),
        "labels": [Path(str(row["assets_dir"])).name for row in summaries],
        "summaries": summaries,
    }
    (outdir / "terminal_climate_replay_summary.json").write_text(
        json.dumps(batch, indent=2, default=_json_default)
    )
    return batch
