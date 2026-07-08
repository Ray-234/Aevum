"""Command-line interface.

    aevum run --preset earthlike --cells 12000 --out out/
    aevum p12 --presets earthlike waterworld --cells 3000 --out out_p12/
    aevum profile-resolution --cells 900 2500 --out out_profile/
    aevum registry [--dump path] [--validate]
    aevum presets
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from aevum.engine import Engine
from aevum.features import build_registry
from aevum.spec.presets import PRESETS, get_preset


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _parse_global_override(raw: str) -> tuple[str, float]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("expected KEY=VALUE")
    key, value = raw.split("=", 1)
    key = key.strip()
    if not key:
        raise argparse.ArgumentTypeError("global override key cannot be empty")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"global override value must be numeric: {value!r}") from exc
    return key, parsed


def cmd_run(args) -> None:
    spec = get_preset(args.preset)
    if args.cells:
        spec.grid_cells = args.cells
    if args.t_end:
        spec.t_end_myr = args.t_end
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"== aevum :: world '{spec.name}' (seed {spec.seed}) ==")
    print(f"   grid {spec.grid_cells} cells, {spec.t_end_myr:.0f} Myr, "
          f"regime {spec.initial_tectonic_regime.value}")
    t0 = time.time()
    eng = Engine.build(spec)
    eng.run(n_frames=args.frames, progress=args.verbose)
    print(f"   deep-time run finished in {time.time() - t0:.1f}s, "
          f"{len(eng.scheduler.history)} macro steps, "
          f"{len(eng.bus.events)} events")

    from aevum.compiler.map_compiler import MapCompiler
    compiler = MapCompiler(eng.world, eng.archive)
    cm = compiler.compile(width=args.hex_w, height=args.hex_h, n_starts=args.starts)
    print(f"   compiled hex map {cm.meta['width']}x{cm.meta['height']}, "
          f"land {cm.meta['land_fraction']:.0%}, fairness {cm.fairness}")

    # validation
    from aevum import validation
    checks = validation.run_all(eng)
    for c in checks:
        print(f"   [{'PASS' if c.passed else 'FAIL'}] {c.name}: {c.detail}")

    # render
    if not args.no_render:
        from aevum import render
        render.render_world(eng.world, outdir)
        render.render_hexmap(cm, outdir)
        render.render_compiler_consistency(cm, outdir)
        render.render_timeline(eng.archive.timeline(), outdir)
        render.render_history(eng.archive, outdir)
        render.render_archive_continuity(eng, outdir)
        print(f"   wrote images to {outdir}/")

    # artefacts
    (outdir / "timeline.json").write_text(
        json.dumps(eng.archive.timeline(), indent=2, default=_json_default))
    # a few example causal explanations (one per start)
    explains = [compiler.explain(cm, r, c) for (r, c) in cm.starts[:3]]
    (outdir / "explain_examples.json").write_text(
        json.dumps(explains, indent=2, default=_json_default))
    (outdir / "lineages.json").write_text(
        json.dumps(eng.archive.lineages(), indent=2, default=_json_default))
    (outdir / "spec.json").write_text(
        json.dumps(spec.to_dict(), indent=2, default=_json_default))

    _print_summary(eng, cm)


def _print_summary(eng, cm) -> None:
    w = eng.world
    print("\n-- world summary --")
    print(f"   mantle T        : {w.g('interior.mantle_temperature'):.0f} K")
    print(f"   tectonic regime : code {w.g('tectonics.regime_code'):.0f}")
    print(f"   CO2             : {w.g('atmosphere.co2')*1e6:.0f} ppm-eq")
    print(f"   O2 fraction     : {w.g('biogeochem.oxygen_fraction'):.4f}")
    print(f"   sea level       : {w.g('ocean.sea_level_m'):.0f} m")
    print(f"   land fraction   : {w.land_fraction():.0%}")
    groups = w.objects.get("biosphere.functional_groups", [])
    alive = [g for g in groups if g['alive']]
    print(f"   functional grps : {len(alive)} alive / {len(groups)} total")
    print(f"   deposits        : {len(w.objects.get('resources.deposits', []))}")
    counts: dict[str, int] = {}
    for e in eng.bus.events:
        counts[e.type] = counts.get(e.type, 0) + 1
    print(f"   events by type  : {counts}")


def cmd_registry(args) -> None:
    reg = build_registry()
    problems = reg.validate()
    loops = reg.feedback_loops()
    print(f"features registered: {len(reg)}")
    print(f"status: {reg.status_summary()}")
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  -", p)
    else:
        print("registry OK (all dependencies resolved)")
    print(f"feedback loops (expected, resolved by iteration): {len(loops)}")
    for lp in loops:
        print("  ~", lp)
    if args.dump:
        reg.dump_yaml(args.dump)
        print(f"dumped to {args.dump}")


def cmd_presets(args) -> None:
    for name, fn in PRESETS.items():
        s = fn()
        print(f"{name:16s} {s.initial_tectonic_regime.value:14s} "
              f"water x{s.composition.water_inventory_earth:<5g} "
              f"land~{s.target_land_fraction:.0%}  {s.notes}")


def cmd_p12(args) -> None:
    from aevum.diagnostics.release_gate import (
        P12RunConfig,
        run_p12_release_gate,
    )

    outdir = Path(args.out)
    config = P12RunConfig(
        presets=tuple(args.presets),
        cells=args.cells,
        t_end_myr=args.t_end,
        frames=args.frames,
        hex_width=args.hex_w,
        hex_height=args.hex_h,
        starts=args.starts,
        render_world_assets=args.render_worlds,
        global_overrides=dict(args.set_global or []),
    )
    t0 = time.time()
    summary = run_p12_release_gate(config, outdir)
    decision = summary["release_decision"]
    print(f"== aevum :: P12 tectonics release gate ==")
    print(f"   presets: {', '.join(config.presets)}")
    print(f"   cells: {config.cells}, frames: {config.frames}, "
          f"runtime {time.time() - t0:.1f}s")
    print(f"   decision: {decision['status']} "
          f"({len(decision['failed_entries'])} failed, "
          f"{len(decision['warned_entries'])} warned)")
    for entry in summary["entries"]:
        gate = entry["release_gate"]
        print(f"   [{gate['status'].upper()}] {entry['preset']}: "
              f"land {entry['land_fraction']:.1%}, "
              f"components {entry['morphology']['land_component_count']}, "
              f"ribbon {entry['morphology']['land_ribbon_fraction_gt_0_5']:.3f}, "
              f"basins {entry['ocean_geography']['basin_count']}")
        for msg in gate["failures"][:3]:
            print(f"      failure: {msg}")
        for msg in gate["warnings"][:3]:
            print(f"      warning: {msg}")
    print(f"   wrote {outdir / 'p12_tectonics_release_summary.json'}")
    print(f"   wrote {outdir / 'p12_preset_matrix_contact_sheet.png'}")


def cmd_profile_resolution(args) -> None:
    from aevum.diagnostics.resolution_profile import (
        ResolutionProfileConfig,
        run_resolution_profile,
    )

    outdir = Path(args.out)
    config = ResolutionProfileConfig(
        preset=args.preset,
        cells=tuple(args.cells),
        t_end_myr=args.t_end,
        frames=args.frames,
        hex_width=args.hex_w,
        hex_height=args.hex_h,
        starts=args.starts,
        compile_map=not args.no_compile,
        tectonic_diagnostics=not args.no_tectonics,
        geomorphology_coverage=args.coverage,
        render_assets=args.render,
        progress=args.verbose,
        projection_cells=tuple(args.project_cells),
    )
    t0 = time.time()
    summary = run_resolution_profile(config, outdir)
    print("== aevum :: resolution profile ==")
    print(f"   preset: {config.preset}")
    print(f"   cells: {', '.join(str(c) for c in config.cells)}")
    print(f"   runtime {time.time() - t0:.1f}s")
    for entry in summary["entries"]:
        stages = entry["stage_seconds"]
        print(
            f"   [{entry['cells']} cells] total {entry['total_seconds']:.2f}s "
            f"run {stages.get('run', 0.0):.2f}s "
            f"compile {stages.get('compile', 0.0):.2f}s "
            f"land {entry['world']['land_fraction']:.1%} "
            f"ribbon {entry['morphology'].get('land_ribbon_fraction_gt_0_5', 0.0):.3f}"
        )
    highres = summary.get("high_resolution", {})
    for proj in highres.get("projection_estimates", [])[:3]:
        minutes = float(proj["estimated_total_seconds"]) / 60.0
        print(f"   projected {proj['target_cells']} cells: ~{minutes:.1f} min")
    print(f"   wrote {outdir / 'resolution_profile_summary.json'}")


def cmd_p107_audit(args) -> None:
    from aevum.diagnostics.p107_audit import (
        P107_PLATE_TERRAIN_MODULES,
        P107AuditConfig,
        P107AuditRun,
        run_p107_audit,
    )

    cells = tuple(int(x) for x in args.cells)
    n_plates = tuple(int(x) for x in args.n_plates)
    if len(n_plates) == 1:
        n_plates = tuple(n_plates[0] for _ in cells)
    if len(n_plates) != len(cells):
        raise SystemExit("--n-plates must have either one value or one value per --cells")
    seeds = tuple(args.seeds or [])
    if seeds and len(seeds) not in {1, len(cells)}:
        raise SystemExit("--seeds must have either one value or one value per --cells")
    if len(seeds) == 1:
        seeds = tuple(seeds[0] for _ in cells)
    runs = tuple(
        P107AuditRun(
            cells=cell_count,
            n_plates=plate_count,
            seed=(None if not seeds else int(seeds[idx])),
            label=f"{cell_count}cells_{plate_count}p",
        )
        for idx, (cell_count, plate_count) in enumerate(zip(cells, n_plates))
    )
    plate_terrain_only = bool(args.plate_terrain_only or args.fast_preview)
    render_world_assets = not args.no_render_world_assets
    render_contact_sheet = not args.no_contact_sheet
    include_earth_reference = not args.no_earth_reference
    if args.fast_preview:
        render_world_assets = False
        render_contact_sheet = False
        include_earth_reference = False
    config = P107AuditConfig(
        preset=args.preset,
        runs=runs,
        t_end_myr=float(args.t_end),
        frames=int(args.frames),
        render_world_assets=render_world_assets,
        render_contact_sheet=render_contact_sheet,
        include_earth_reference=include_earth_reference,
        enable_ranked_plate_policy=not args.disable_ranked_plate_policy,
        enable_boundary_province_response=not args.disable_boundary_province_response,
        enable_p108_boundary_width_guard=not args.disable_p108_boundary_width_guard,
        enable_p108_high_mountain_coherence=not args.disable_p108_high_mountain_coherence,
        enabled_modules=(
            P107_PLATE_TERRAIN_MODULES if plate_terrain_only else None
        ),
        global_overrides=dict(args.set_global or []),
    )
    outdir = Path(args.out)
    t0 = time.time()
    summary = run_p107_audit(config, outdir)
    print("== aevum :: P107 terminal audit ladder ==")
    print(f"   preset: {config.preset}")
    print(
        "   modules: "
        + (
            "full"
            if config.enabled_modules is None
            else ", ".join(config.enabled_modules)
        )
    )
    if args.fast_preview:
        print("   mode: fast preview (no resources; no render/reference)")
    print(f"   runtime {time.time() - t0:.1f}s")
    for entry in summary["entries"]:
        metrics = entry["metrics"]
        stages = entry.get("stage_seconds", {})
        audit_seconds = float(stages.get("terminal_audit_write", 0.0))
        print(
            f"   [{entry['label']}] cells {entry['cells']} plates {entry['n_plates']} "
            f"run {entry['run_seconds']:.1f}s audit {audit_seconds:.1f}s "
            f"land {metrics['land_fraction']:.1%} "
            f"active/major/minor/micro "
            f"{metrics['terminal_active_plate_count']}/"
            f"{metrics['terminal_major_plate_count']}/"
            f"{metrics['terminal_minor_plate_count']}/"
            f"{metrics['terminal_microplate_count']}"
        )
        top_modules = entry.get("module_seconds", {}).get("top", [])[:3]
        if top_modules:
            module_text = ", ".join(
                f"{row['module']}={float(row['seconds']):.1f}s"
                for row in top_modules
            )
            print(f"      module time: {module_text}")
        terrain_profile = (
            entry.get("terrain_internal_profile")
            or metrics.get("terrain_internal_profile", {})
        )
        terrain_top = terrain_profile.get("top", []) if terrain_profile else []
        if terrain_top:
            terrain_text = ", ".join(
                f"{row['stage']}={float(row['seconds']):.2f}s"
                for row in terrain_top[:3]
            )
            print(f"      terrain profile: {terrain_text}")
        semantic_profile = (
            terrain_profile.get("subprofiles", {})
            .get("semantic_object_build", {})
            if terrain_profile else {}
        )
        semantic_top = semantic_profile.get("top", []) if semantic_profile else []
        if semantic_top:
            semantic_text = ", ".join(
                f"{row['stage']}={float(row['seconds']):.2f}s"
                for row in semantic_top[:3]
            )
            print(f"      semantic profile: {semantic_text}")
        terminal_bathymetry_profile = (
            terrain_profile.get("subprofiles", {})
            .get("terminal_bathymetry_polish", {})
            if terrain_profile else {}
        )
        terminal_bathymetry_top = (
            terminal_bathymetry_profile.get("top", [])
            if terminal_bathymetry_profile else []
        )
        if terminal_bathymetry_top:
            terminal_bathymetry_text = ", ".join(
                f"{row['stage']}={float(row['seconds']):.2f}s"
                for row in terminal_bathymetry_top[:3]
            )
            print(f"      bathymetry profile: {terminal_bathymetry_text}")
        print(f"      wrote {entry['outdir']}")
    print(f"   wrote {outdir / 'p107_audit_summary.json'}")


def cmd_p110b_seed_sweep(args) -> None:
    from aevum.diagnostics.p110b_seed_sweep import (
        P110BSeedSweepThresholds,
        summarize_p110b_seed_sweep,
        write_p110b_seed_sweep_summary,
    )

    thresholds = P110BSeedSweepThresholds(
        min_sample_size=int(args.min_sample_size),
        max_soft_warning_rate=float(args.max_soft_warning_rate),
        max_median_largest_land_share=float(args.max_median_largest_share),
        max_p90_largest_land_share=float(args.max_p90_largest_share),
    )
    summary = summarize_p110b_seed_sweep(args.inputs, thresholds=thresholds)
    out_path = write_p110b_seed_sweep_summary(summary, args.out)
    aggregate = summary["aggregate"]
    largest = aggregate["largest_land_component_share"]
    rates = aggregate["rates"]
    print("== aevum :: P110B seed-sweep diagnostic ==")
    print(f"   inputs: {len(summary['input_paths'])}, runs: {summary['run_count']}")
    print(
        f"   P109/P110A threshold pass: {aggregate['threshold_pass_count']}/"
        f"{summary['run_count']}"
    )
    print(
        f"   P110B visual candidates: "
        f"{aggregate['p110b_visual_candidate_count']}/{summary['run_count']}"
    )
    print(
        f"   P111 modern-planform candidates: "
        f"{aggregate.get('p111_modern_planform_candidate_count', 0)}/"
        f"{summary['run_count']}"
    )
    print(
        f"   largest land share median/p90/max: "
        f"{(largest.get('median') or 0.0):.3f}/"
        f"{(largest.get('p90') or 0.0):.3f}/"
        f"{(largest.get('max') or 0.0):.3f}"
    )
    print(
        f"   largest soft-warning rate: "
        f"{rates['largest_land_component_soft_warning_rate']:.1%}"
    )
    flags = aggregate["distribution_flags"]
    if flags:
        print(f"   distribution flags: {', '.join(flags)}")
    else:
        print("   distribution flags: none")
    p111_flags = aggregate.get("p111_distribution_flags", [])
    if p111_flags:
        print(f"   P111 distribution flags: {', '.join(p111_flags)}")
    else:
        print("   P111 distribution flags: none")
    print(f"   wrote {out_path}")


def cmd_p107_render_arrays(args) -> None:
    from aevum.diagnostics.p107_array_render import render_p107_array_assets

    summary = render_p107_array_assets(
        args.input,
        args.out,
        width=int(args.width),
        height=int(args.height),
    )
    assets = summary["assets"]
    p110a = summary.get("p110a_summary", {})
    warnings = summary.get("p110a_warning_flags", [])
    print("== aevum :: P107 array renderer ==")
    print(f"   source: {summary['source_metrics']}")
    print(f"   cells: {summary['cells']}, assets: {len(assets)}")
    print(
        f"   P110A largest/second/third: "
        f"{p110a.get('largest_land_component_share', 0.0):.3f}/"
        f"{p110a.get('second_land_component_share', 0.0):.3f}/"
        f"{p110a.get('third_land_component_share', 0.0):.3f}"
    )
    if warnings:
        print(f"   P110A warnings: {', '.join(str(x) for x in warnings)}")
    print(f"   wrote {Path(args.out) / 'p107_array_render_summary.json'}")


def cmd_selected_snapshot_refine(args) -> None:
    from aevum.diagnostics.selected_snapshot_refinement import (
        SelectedSnapshotRefinementConfig,
        refine_selected_snapshot,
    )

    summary = refine_selected_snapshot(
        SelectedSnapshotRefinementConfig(
            source=args.input,
            outdir=args.out,
            target_cells=int(args.target_cells),
            width=int(args.width),
            height=int(args.height),
            interpolation_k=int(args.interpolation_k),
            detail_seed=int(args.detail_seed),
            detail_strength=float(args.detail_strength),
            allow_process_islands=bool(args.allow_process_islands),
            render_groups=tuple(args.render_groups or ("all",)),
        )
    )
    print("== aevum :: selected snapshot refinement ==")
    print(f"   source: {summary['source_metrics']}")
    print(
        f"   cells: {summary['source_cells']} -> {summary['target_cells']}, "
        f"land {summary['land_fraction_parent']:.3f} -> "
        f"{summary['land_fraction_refined']:.3f}"
    )
    print(
        f"   sign flips: {summary['land_ocean_sign_flip_fraction']:.6f}, "
        f"delta p95 land/ocean: "
        f"{summary['detail_delta_land_abs_p95_m']:.1f}/"
        f"{summary['detail_delta_ocean_abs_p95_m']:.1f} m"
    )
    if summary.get("allow_process_islands"):
        print(
            f"   process islands: "
            f"{summary.get('process_island_promoted_cell_fraction_parent_ocean', 0.0):.6f} "
            f"of parent ocean, atoll islets "
            f"{summary.get('atoll_islet_promoted_cell_fraction_parent_ocean', 0.0):.6f}"
        )
    print(f"   render groups: {', '.join(summary.get('render_groups', ['all']))}")
    print(f"   wrote {Path(args.out) / 'selected_snapshot_refinement_metrics.json'}")


def cmd_selected_snapshot_render_groups(args) -> None:
    from aevum.diagnostics.selected_snapshot_refinement import (
        render_selected_snapshot_refinement_assets,
    )

    summary = render_selected_snapshot_refinement_assets(
        args.input,
        render_groups=tuple(args.render_groups or ("all",)),
        width=args.width,
        height=args.height,
        outdir=args.out,
    )
    assets = summary["assets"]["refinement"]
    print("== aevum :: selected snapshot render groups ==")
    print(f"   source: {summary['source_metrics']}")
    print(f"   render groups: {', '.join(summary['render_groups'])}")
    print(f"   refinement assets: {len(assets)}")
    rendered_dir = Path(args.out) if args.out else Path(summary["source_metrics"]).parent / "rendered"
    print(f"   wrote {rendered_dir}")


def cmd_p107_compare(args) -> None:
    from aevum.diagnostics.p107_equivalence import (
        DEFAULT_METRIC_SKIP_KEYS,
        compare_p107_outputs,
        write_p107_equivalence_report,
    )

    skip_keys = (
        tuple(args.metric_skip_key)
        if args.metric_skip_key
        else DEFAULT_METRIC_SKIP_KEYS
    )
    report = compare_p107_outputs(
        args.baseline,
        args.candidate,
        skip_top_keys=skip_keys,
        float_atol=float(args.float_atol),
    )
    out_path = write_p107_equivalence_report(report, args.out)
    arrays = report["arrays"]
    metrics = report["metrics"]
    status = "PASS" if report["equivalent"] else "FAIL"
    print("== aevum :: P107 output equivalence ==")
    print(f"   status: {status}")
    print(f"   baseline: {report['baseline']['run_dir']}")
    print(f"   candidate: {report['candidate']['run_dir']}")
    print(
        f"   arrays: common {arrays['common_count']}, "
        f"changed {arrays['changed_count']}, "
        f"missing {len(arrays['missing_keys'])}, "
        f"extra {len(arrays['extra_keys'])}, "
        f"dtype mismatches {len(arrays['dtype_mismatches'])}"
    )
    print(
        f"   metrics: changed {metrics['changed_count']}, "
        f"max float diff {metrics['max_float_diff']['diff']} "
        f"at {metrics['max_float_diff']['path']!r}"
    )
    print(f"   wrote {out_path}")
    if not report["equivalent"]:
        raise SystemExit(1)


def _parse_terminal_climate_job(raw: str):
    parts = raw.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected PRESET:LABEL:SEED")
    preset, label, seed_raw = parts
    if preset not in PRESETS:
        raise argparse.ArgumentTypeError(
            f"unknown preset {preset!r}; options: {sorted(PRESETS)}")
    try:
        seed = int(seed_raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("SEED must be an integer") from exc
    from aevum.diagnostics.terminal_climate_biome import TerminalClimateJob
    return TerminalClimateJob(preset, label, seed)


def cmd_terminal_climate_biome(args) -> None:
    from aevum.diagnostics.terminal_climate_biome import (
        DEFAULT_SIX_WORLD_JOBS,
        TerminalClimateConfig,
        run_terminal_climate_biome_batch,
    )

    jobs = tuple(args.job) if args.job else DEFAULT_SIX_WORLD_JOBS
    config = TerminalClimateConfig(
        jobs=jobs,
        cells=int(args.cells),
        t_end_myr=float(args.t_end),
        frames=int(args.frames),
        max_workers=int(args.max_workers),
        render_assets=not args.no_render,
    )
    outdir = Path(args.out)
    t0 = time.time()
    summary = run_terminal_climate_biome_batch(config, outdir)
    print("== aevum :: terminal climate/biome batch ==")
    print(f"   jobs: {len(jobs)}, cells: {config.cells}, "
          f"workers: {config.max_workers}, runtime {time.time() - t0:.1f}s")
    for entry in summary["summaries"]:
        label = Path(entry["assets_dir"]).name
        print(
            f"   [{label}] land {entry['land_fraction']:.1%} "
            f"T {entry['mean_temperature_C']:.1f} C "
            f"P {entry['mean_precipitation_mm_yr']:.0f} mm/yr"
        )
    print(f"   wrote {outdir / 'terminal_climate_biome_summary.json'}")


def cmd_earth_climate_reference(args) -> None:
    from aevum.diagnostics.earth_climate_reference import (
        EarthClimateReferenceConfig,
        run_earth_climate_reference,
    )

    manifest_out = Path(args.manifest_out) if args.manifest_out else None
    config = EarthClimateReferenceConfig(
        cells=tuple(int(x) for x in args.cells),
        width=int(args.width),
        height=int(args.height),
        render_assets=not args.no_render,
        manifest_out=manifest_out,
        include_worldclim=not args.no_worldclim,
        include_koppen=not args.no_koppen,
        include_noaa_psl=not args.no_noaa_psl,
        include_sst=not args.no_sst,
        include_ocean_currents=not args.no_ocean_currents,
        include_etopo2022=not args.no_etopo2022,
        include_land_cover=not args.no_land_cover,
        download_oscar=bool(args.download_oscar),
        download_land_cover=bool(args.download_land_cover),
        climatology_start_year=int(args.climatology_start_year),
        climatology_end_year=int(args.climatology_end_year),
    )
    outdir = Path(args.out)
    t0 = time.time()
    summary = run_earth_climate_reference(config, outdir)
    print("== aevum :: Earth climate reference ==")
    print(
        f"   cells: {', '.join(str(c) for c in config.cells)}, "
        f"runtime {time.time() - t0:.1f}s"
    )
    for entry in summary["entries"]:
        metrics = entry["metrics"]["elevation"]
        print(
            f"   [{entry['cells']} cells] land {metrics['land_fraction']:.1%} "
            f"land mean {metrics['land_elevation_mean_m']:.0f} m "
            f"ocean mean depth {metrics['ocean_depth_mean_m']:.0f} m"
        )
        if "worldclim" in entry["metrics"]:
            climate = entry["metrics"]["worldclim"]
            print(
                f"      WorldClim land T "
                f"{climate['land_annual_temperature_mean_C']:.1f} C, "
                f"P {climate['land_annual_precip_mean_mm']:.0f} mm/yr"
            )
        if "noaa_oisst_v2" in entry["metrics"]:
            sst = entry["metrics"]["noaa_oisst_v2"]
            print(
                f"      NOAA OISST SST mean "
                f"{sst['annual_sst_mean_C']:.1f} C, "
                f"tropics {sst['tropical_sst_mean_C']:.1f} C"
            )
        if "noaa_psl_ncep" in entry["metrics"]:
            ncep = entry["metrics"]["noaa_psl_ncep"]
            print(
                f"      NOAA PSL wind p90 "
                f"{ncep['seasonal_wind_speed_p90_m_s']:.1f} m/s, "
                f"SLP anomaly p90 "
                f"{ncep['seasonal_slp_anomaly_abs_p90_hPa']:.1f} hPa"
            )
        if "noaa_aoml_drifter_current_v3" in entry["metrics"]:
            currents = entry["metrics"]["noaa_aoml_drifter_current_v3"]
            print(
                f"      NOAA/AOML current speed p90 "
                f"{currents['current_speed_p90_m_s']:.2f} m/s"
            )
        if "nasa_jpl_oscar_monthly" in entry["metrics"]:
            currents = entry["metrics"]["nasa_jpl_oscar_monthly"]
            print(
                f"      NASA/JPL OSCAR annual current speed p90 "
                f"{currents['annual_current_speed_p90_m_s']:.2f} m/s, "
                f"monthly p90 {currents['monthly_current_speed_p90_m_s']:.2f} m/s"
            )
        if "esa_cci_land_cover" in entry["metrics"]:
            land_cover = entry["metrics"]["esa_cci_land_cover"]
            print(
                f"      ESA CCI LC forest "
                f"{land_cover['forest_area_fraction']:.1%}, "
                f"crop {land_cover['cropland_area_fraction']:.1%}, "
                f"water {land_cover['water_area_fraction']:.1%}"
            )
    print(f"   wrote {outdir / 'earth_climate_reference_summary.json'}")
    if manifest_out:
        print(f"   wrote {manifest_out}")


def cmd_earth_climate_compare(args) -> None:
    from aevum.diagnostics.earth_climate_comparison import (
        EarthClimateComparisonConfig,
        run_earth_climate_comparison,
    )

    config = EarthClimateComparisonConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
        render_contact_sheet=not args.no_contact_sheet,
    )
    summary = run_earth_climate_comparison(config)
    print("== aevum :: Earth climate comparison ==")
    print(f"   runs: {summary['run_count']}")
    print(f"   earthlike flagged: {summary['earthlike_flagged_count']}")
    for entry in summary["entries"]:
        print(
            f"   [{entry['mode']}] {entry['label']}: "
            f"score {entry['earth_distance_score']:.2f}, "
            f"flags {len(entry['flags'])}"
        )
        for flag in entry["flags"][:4]:
            print(f"      {flag}")
    print(f"   wrote {Path(args.out) / 'earth_climate_comparison_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_comparison_metrics.csv'}")
    if summary.get("contact_sheet"):
        print(f"   wrote {summary['contact_sheet']}")


def cmd_earth_climate_fit_report(args) -> None:
    from aevum.diagnostics.earth_climate_fitting import (
        EarthClimateFittingConfig,
        run_earth_climate_fitting_report,
    )

    config = EarthClimateFittingConfig(
        comparison_summary_json=Path(args.comparison_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_fitting_report(config)
    print("== aevum :: Earth climate fitting report ==")
    print(f"   earthlike runs: {report['earthlike_run_count']}")
    print(f"   conclusion: {report['dominant_conclusion']}")
    guard = report.get("guardrail_assessment", {})
    print(
        f"   guardrail verdict: {report.get('overall_verdict', 'unknown')} "
        f"({guard.get('failure_count', 0)} failures, "
        f"{guard.get('warning_count', 0)} warnings)"
    )
    for phase, row in report["phase_assessment"].items():
        print(
            f"   [{row['priority']}] {phase}: "
            f"{row['status']} score {row['score']:.2f}"
        )
    for warning in guard.get("warnings", [])[:4]:
        print(f"      warning: {warning['label']} {warning['metric']}")
    for failure in guard.get("failures", [])[:4]:
        print(f"      failure: {failure['label']} {failure['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_fitting_report.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_fitting_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_fitting_runs.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_fitting_levers.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_guardrails.csv'}")
    if args.fail_on_guardrail and report.get("overall_verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_pattern_gate(args) -> None:
    from aevum.diagnostics.earth_climate_pattern_gate import (
        EarthClimatePatternGateConfig,
        run_earth_climate_pattern_gate,
    )

    config = EarthClimatePatternGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_pattern_gate(config)
    print("== aevum :: Earth climate pattern gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_pattern_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_pattern_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_pattern_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_pattern_checks.csv'}")
    if args.fail_on_pattern and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_biome_gate(args) -> None:
    from aevum.diagnostics.earth_climate_biome_gate import (
        EarthClimateBiomeGateConfig,
        run_earth_climate_biome_gate,
    )

    config = EarthClimateBiomeGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_biome_gate(config)
    print("== aevum :: Earth climate biome gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_biome_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_biome_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_biome_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_biome_checks.csv'}")
    if args.fail_on_biome and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_spatial_biome_gate(args) -> None:
    from aevum.diagnostics.earth_climate_spatial_biome_gate import (
        EarthClimateSpatialBiomeGateConfig,
        run_earth_climate_spatial_biome_gate,
    )

    config = EarthClimateSpatialBiomeGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_spatial_biome_gate(config)
    print("== aevum :: Earth climate spatial biome gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_spatial_biome_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_spatial_biome_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_spatial_biome_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_spatial_biome_checks.csv'}")
    if args.fail_on_spatial_biome and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_seasonal_subtype_gate(args) -> None:
    from aevum.diagnostics.earth_climate_seasonal_subtype_gate import (
        EarthClimateSeasonalSubtypeGateConfig,
        run_earth_climate_seasonal_subtype_gate,
    )

    config = EarthClimateSeasonalSubtypeGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_seasonal_subtype_gate(config)
    print("== aevum :: Earth climate seasonal subtype gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_subtype_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_subtype_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_subtype_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_subtype_checks.csv'}")
    if args.fail_on_seasonal_subtype and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_mountain_zonation_gate(args) -> None:
    from aevum.diagnostics.earth_climate_mountain_zonation_gate import (
        EarthClimateMountainZonationGateConfig,
        run_earth_climate_mountain_zonation_gate,
    )

    config = EarthClimateMountainZonationGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_mountain_zonation_gate(config)
    print("== aevum :: Earth climate mountain zonation gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_mountain_zonation_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_mountain_zonation_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_mountain_zonation_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_mountain_zonation_checks.csv'}")
    if args.fail_on_mountain_zonation and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_windward_leeward_gate(args) -> None:
    from aevum.diagnostics.earth_climate_windward_leeward_gate import (
        EarthClimateWindwardLeewardGateConfig,
        run_earth_climate_windward_leeward_gate,
    )

    config = EarthClimateWindwardLeewardGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_windward_leeward_gate(config)
    print("== aevum :: Earth climate windward/leeward gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_windward_leeward_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_windward_leeward_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_windward_leeward_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_windward_leeward_checks.csv'}")
    if args.fail_on_windward_leeward and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_monsoon_moisture_gate(args) -> None:
    from aevum.diagnostics.earth_climate_monsoon_moisture_gate import (
        EarthClimateMonsoonMoistureGateConfig,
        run_earth_climate_monsoon_moisture_gate,
    )

    config = EarthClimateMonsoonMoistureGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_monsoon_moisture_gate(config)
    print("== aevum :: Earth climate monsoon/moisture gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_monsoon_moisture_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_monsoon_moisture_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_monsoon_moisture_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_monsoon_moisture_checks.csv'}")
    if args.fail_on_monsoon_moisture and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_hydro_region_gate(args) -> None:
    from aevum.diagnostics.earth_climate_hydro_region_gate import (
        EarthClimateHydroRegionGateConfig,
        run_earth_climate_hydro_region_gate,
    )

    config = EarthClimateHydroRegionGateConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
        render_contact_sheets=not bool(getattr(args, "no_contact_sheet", False)),
    )
    report = run_earth_climate_hydro_region_gate(config)
    print("== aevum :: Earth climate hydroclimate-region gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_hydro_region_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_hydro_region_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_hydro_region_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_hydro_region_checks.csv'}")
    if report.get("contact_sheet_count", 0):
        print(
            f"   wrote {report['contact_sheet_count']} contact sheets under "
            f"{Path(args.out) / 'contact_sheets'}"
        )
    if args.fail_on_hydro_regions and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_moisture_flow_gate(args) -> None:
    from aevum.diagnostics.earth_climate_moisture_flow_gate import (
        EarthClimateMoistureFlowGateConfig,
        run_earth_climate_moisture_flow_gate,
    )

    config = EarthClimateMoistureFlowGateConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
        render_contact_sheets=not bool(getattr(args, "no_contact_sheet", False)),
    )
    report = run_earth_climate_moisture_flow_gate(config)
    print("== aevum :: Earth climate moisture-flow-network gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_flow_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_flow_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_flow_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_flow_checks.csv'}")
    if report.get("contact_sheet_count", 0):
        print(
            f"   wrote {report['contact_sheet_count']} contact sheets under "
            f"{Path(args.out) / 'contact_sheets'}"
        )
    if args.fail_on_moisture_flow and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_moisture_response_gate(args) -> None:
    from aevum.diagnostics.earth_climate_moisture_response_gate import (
        EarthClimateMoistureResponseGateConfig,
        run_earth_climate_moisture_response_gate,
    )

    config = EarthClimateMoistureResponseGateConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
        render_contact_sheets=not bool(getattr(args, "no_contact_sheet", False)),
    )
    report = run_earth_climate_moisture_response_gate(config)
    print("== aevum :: Earth climate moisture-flow precipitation-response gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_response_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_response_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_response_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_moisture_response_checks.csv'}")
    if report.get("contact_sheet_count", 0):
        print(
            f"   wrote {report['contact_sheet_count']} contact sheets under "
            f"{Path(args.out) / 'contact_sheets'}"
        )
    if args.fail_on_moisture_response and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_receiver_catchment_gate(args) -> None:
    from aevum.diagnostics.earth_climate_receiver_catchment_gate import (
        EarthClimateReceiverCatchmentGateConfig,
        run_earth_climate_receiver_catchment_gate,
    )

    config = EarthClimateReceiverCatchmentGateConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_receiver_catchment_gate(config)
    print("== aevum :: Earth climate receiver-catchment gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(
        f"   wrote {Path(args.out) / 'earth_climate_receiver_catchment_gate_summary.json'}"
    )
    print(
        f"   wrote {Path(args.out) / 'earth_climate_receiver_catchment_gate_report.md'}"
    )
    print(f"   wrote {Path(args.out) / 'earth_climate_receiver_catchment_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_receiver_catchment_checks.csv'}")
    if args.fail_on_receiver_catchment and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_circulation_layout_gate(args) -> None:
    from aevum.diagnostics.earth_climate_circulation_layout_gate import (
        EarthClimateCirculationLayoutGateConfig,
        run_earth_climate_circulation_layout_gate,
    )

    config = EarthClimateCirculationLayoutGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_circulation_layout_gate(config)
    print("== aevum :: Earth climate circulation-layout gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(
        f"   wrote {Path(args.out) / 'earth_climate_circulation_layout_gate_summary.json'}"
    )
    print(
        f"   wrote {Path(args.out) / 'earth_climate_circulation_layout_gate_report.md'}"
    )
    print(f"   wrote {Path(args.out) / 'earth_climate_circulation_layout_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_circulation_layout_checks.csv'}")
    if args.fail_on_circulation_layout and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_real_earth_wind_replay(args) -> None:
    from aevum.diagnostics.real_earth_wind_replay import (
        RealEarthWindReplayConfig,
        run_real_earth_wind_replay,
    )

    config = RealEarthWindReplayConfig(
        earth_reference_npz=Path(args.earth_reference),
        replay_arrays_npz=Path(args.replay_arrays),
        outdir=Path(args.out),
    )
    summary = run_real_earth_wind_replay(config)
    metrics = summary["metrics"]
    print("== aevum :: real-Earth wind replay ==")
    print(
        "   seasonal wind speed MAE: "
        f"{metrics['seasonal_speed_mae_m_s']:.2f} m/s; "
        "direction cosine p50: "
        f"{metrics['direction_cosine_p50']:.2f}"
    )
    print(f"   wrote {Path(args.out) / 'real_earth_wind_replay_summary.json'}")
    print(f"   wrote {Path(args.out) / 'real_earth_wind_replay_contact_sheet.png'}")


def cmd_earth_climate_ocean_spatial_gate(args) -> None:
    from aevum.diagnostics.earth_climate_ocean_spatial_gate import (
        EarthClimateOceanSpatialGateConfig,
        run_earth_climate_ocean_spatial_gate,
    )

    config = EarthClimateOceanSpatialGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_ocean_spatial_gate(config)
    print("== aevum :: Earth climate ocean-current/SST spatial gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(f"   wrote {Path(args.out) / 'earth_climate_ocean_spatial_gate_summary.json'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_ocean_spatial_gate_report.md'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_ocean_spatial_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_ocean_spatial_checks.csv'}")
    if args.fail_on_ocean_spatial and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_coupled_consistency_gate(args) -> None:
    from aevum.diagnostics.earth_climate_coupled_consistency_gate import (
        EarthClimateCoupledConsistencyGateConfig,
        run_earth_climate_coupled_consistency_gate,
    )

    config = EarthClimateCoupledConsistencyGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_coupled_consistency_gate(config)
    print("== aevum :: Earth climate coupled-consistency gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(
        f"   wrote {Path(args.out) / 'earth_climate_coupled_consistency_gate_summary.json'}"
    )
    print(
        f"   wrote {Path(args.out) / 'earth_climate_coupled_consistency_gate_report.md'}"
    )
    print(f"   wrote {Path(args.out) / 'earth_climate_coupled_consistency_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_coupled_consistency_checks.csv'}")
    if args.fail_on_coupled_consistency and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_seasonal_hydro_placement_gate(args) -> None:
    from aevum.diagnostics.earth_climate_seasonal_hydro_placement_gate import (
        EarthClimateSeasonalHydroPlacementGateConfig,
        run_earth_climate_seasonal_hydro_placement_gate,
    )

    config = EarthClimateSeasonalHydroPlacementGateConfig(
        earth_reference_npz=Path(args.earth_reference),
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_seasonal_hydro_placement_gate(config)
    print("== aevum :: Earth climate seasonal-hydro placement gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(
        f"   wrote {Path(args.out) / 'earth_climate_seasonal_hydro_placement_gate_summary.json'}"
    )
    print(
        f"   wrote {Path(args.out) / 'earth_climate_seasonal_hydro_placement_gate_report.md'}"
    )
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_hydro_placement_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_seasonal_hydro_placement_checks.csv'}")
    if args.fail_on_seasonal_hydro_placement and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_earth_climate_coupling_convergence_gate(args) -> None:
    from aevum.diagnostics.earth_climate_coupling_convergence_gate import (
        EarthClimateCouplingConvergenceGateConfig,
        run_earth_climate_coupling_convergence_gate,
    )

    config = EarthClimateCouplingConvergenceGateConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
    )
    report = run_earth_climate_coupling_convergence_gate(config)
    print("== aevum :: Earth climate coupling-convergence gate ==")
    print(
        f"   verdict: {report['verdict']} "
        f"({report['failure_count']} failures, "
        f"{report['warning_count']} warnings, "
        f"{report['skipped_count']} skipped)"
    )
    for failure in report.get("failures", [])[:6]:
        print(f"      failure: {failure['label']} {failure['group']} {failure['metric']}")
    for warning in report.get("warnings", [])[:6]:
        print(f"      warning: {warning['label']} {warning['group']} {warning['metric']}")
    print(
        f"   wrote {Path(args.out) / 'earth_climate_coupling_convergence_gate_summary.json'}"
    )
    print(
        f"   wrote {Path(args.out) / 'earth_climate_coupling_convergence_gate_report.md'}"
    )
    print(f"   wrote {Path(args.out) / 'earth_climate_coupling_convergence_metrics.csv'}")
    print(f"   wrote {Path(args.out) / 'earth_climate_coupling_convergence_checks.csv'}")
    if args.fail_on_coupling_convergence and report.get("verdict") == "fail":
        raise SystemExit(2)


def cmd_terminal_climate_replay(args) -> None:
    from aevum.diagnostics.terminal_climate_replay import (
        TerminalClimateReplayConfig,
        run_terminal_climate_replay,
    )

    config = TerminalClimateReplayConfig(
        terminal_summary_json=Path(args.terminal_summary),
        outdir=Path(args.out),
        labels=tuple(args.label or ()),
        render_assets=not args.no_render,
    )
    summary = run_terminal_climate_replay(config)
    print("== aevum :: terminal climate replay ==")
    print(f"   jobs: {summary['job_count']}")
    for entry in summary["summaries"]:
        print(
            f"   {Path(entry['assets_dir']).name}: "
            f"T {entry['mean_temperature_C']:.2f} C, "
            f"land precip p50 {entry['land_precip_p50_mm_yr']}"
        )
    print(f"   wrote {Path(args.out) / 'terminal_climate_replay_summary.json'}")


def cmd_real_earth_climate_replay(args) -> None:
    from aevum.diagnostics.real_earth_climate_replay import (
        RealEarthClimateReplayConfig,
        run_real_earth_climate_replay,
    )

    config = RealEarthClimateReplayConfig(
        earth_reference_npz=Path(args.earth_reference),
        outdir=Path(args.out),
        preset=args.preset,
        seed=int(args.seed),
        render_assets=not args.no_render,
    )
    summary = run_real_earth_climate_replay(config)
    residuals = summary["residuals"]
    dynamics = summary["dynamics_residuals"]
    print("== aevum :: real-Earth climate replay ==")
    print(f"   cells: {summary['cells']}")
    print(
        "   T MAE "
        f"{residuals['surface_temperature_mae_C']:.2f} C, "
        f"land precip MAE {residuals['land_precip_mae_mm_yr']:.1f} mm/yr"
    )
    print(
        "   wind p90 replay/earth "
        f"{dynamics['annual_wind_speed_replay_p90_m_s']:.2f}/"
        f"{dynamics['annual_wind_speed_earth_p90_m_s']:.2f} m/s"
    )
    print(f"   validation: {'PASS' if summary['validation']['passed'] else 'FAIL'}")
    print(f"   wrote {Path(args.out) / 'real_earth_climate_replay_summary.json'}")


def main() -> None:
    ap = argparse.ArgumentParser(prog="aevum",
                                 description="planetary deep-time evolution engine")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run a world end to end and compile a map")
    r.add_argument("--preset", default="earthlike", choices=sorted(PRESETS))
    r.add_argument("--cells", type=int, default=8000)
    r.add_argument("--t-end", type=float, default=None, dest="t_end")
    r.add_argument("--frames", type=int, default=18)
    r.add_argument("--hex-w", type=int, default=96)
    r.add_argument("--hex-h", type=int, default=48)
    r.add_argument("--starts", type=int, default=6)
    r.add_argument("--out", default="out")
    r.add_argument("--no-render", action="store_true")
    r.add_argument("-v", "--verbose", action="store_true")
    r.set_defaults(func=cmd_run)

    g = sub.add_parser("registry", help="inspect / dump the feature registry")
    g.add_argument("--dump", default=None)
    g.add_argument("--validate", action="store_true")
    g.set_defaults(func=cmd_registry)

    p = sub.add_parser("presets", help="list baseline worlds")
    p.set_defaults(func=cmd_presets)

    p12 = sub.add_parser("p12", help="run the tectonics release-gate matrix")
    p12.add_argument("--presets", nargs="+", default=[
        "earthlike", "waterworld", "arid", "stagnant_lid", "tidally_locked", "frozen",
    ], choices=sorted(PRESETS))
    p12.add_argument("--cells", type=int, default=3000)
    p12.add_argument("--t-end", type=float, default=None, dest="t_end")
    p12.add_argument("--frames", type=int, default=4)
    p12.add_argument("--hex-w", type=int, default=64)
    p12.add_argument("--hex-h", type=int, default=32)
    p12.add_argument("--starts", type=int, default=4)
    p12.add_argument("--out", default="out_p12")
    p12.add_argument(
        "--set-global",
        action="append",
        type=_parse_global_override,
        default=[],
        metavar="KEY=VALUE",
        help="set a numeric world global before running each P12 world",
    )
    p12.add_argument("--render-worlds", action="store_true",
                     help="also render full per-preset world assets")
    p12.set_defaults(func=cmd_p12)

    prof = sub.add_parser(
        "profile-resolution",
        help="profile a resolution ladder before high-resolution geomorphology runs",
    )
    prof.add_argument("--preset", default="earthlike", choices=sorted(PRESETS))
    prof.add_argument("--cells", nargs="+", type=int, default=[900, 2500])
    prof.add_argument("--t-end", type=float, default=None, dest="t_end")
    prof.add_argument("--frames", type=int, default=4)
    prof.add_argument("--hex-w", type=int, default=64)
    prof.add_argument("--hex-h", type=int, default=32)
    prof.add_argument("--starts", type=int, default=4)
    prof.add_argument("--out", default="out_resolution_profile")
    prof.add_argument("--no-compile", action="store_true",
                      help="skip hex compilation for large scaling probes")
    prof.add_argument("--no-tectonics", action="store_true",
                      help="skip tectonic diagnostics and morphology metrics")
    prof.add_argument("--coverage", action="store_true",
                      help="also compute integrated geomorphology coverage metrics")
    prof.add_argument("--render", action="store_true",
                      help="render per-resolution world assets")
    prof.add_argument("--project-cells", nargs="+", type=int,
                      default=[8000, 24000, 72000],
                      help="cell counts to extrapolate for high-resolution audits")
    prof.add_argument("-v", "--verbose", action="store_true")
    prof.set_defaults(func=cmd_profile_resolution)

    p107 = sub.add_parser(
        "p107-audit",
        help="run the P107 terminal plate/boundary high-resolution audit ladder",
    )
    p107.add_argument("--preset", default="earthlike", choices=sorted(PRESETS))
    p107.add_argument("--cells", nargs="+", type=int, default=[8000, 24000])
    p107.add_argument("--n-plates", nargs="+", type=int, default=[36, 60],
                      dest="n_plates")
    p107.add_argument("--seeds", nargs="+", type=int, default=[])
    p107.add_argument("--t-end", type=float, default=4500.0, dest="t_end")
    p107.add_argument("--frames", type=int, default=5)
    p107.add_argument("--out", default="out_p107_audit")
    p107.add_argument(
        "--set-global",
        action="append",
        type=_parse_global_override,
        default=[],
        metavar="KEY=VALUE",
        help="set a numeric world global before running each P107 world",
    )
    p107.add_argument("--no-render-world-assets", action="store_true")
    p107.add_argument("--no-contact-sheet", action="store_true")
    p107.add_argument("--no-earth-reference", action="store_true")
    p107.add_argument(
        "--plate-terrain-only",
        action="store_true",
        help=(
            "skip resource genesis while retaining climate/biogeochem/biosphere "
            "feedbacks needed by terrain; this is not a full release audit"
        ),
    )
    p107.add_argument(
        "--fast-preview",
        action="store_true",
        help=(
            "shortcut for plate-terrain-only plus no render assets, no contact "
            "sheet, and no Earth-reference extraction"
        ),
    )
    p107.add_argument("--disable-ranked-plate-policy", action="store_true",
                      help="run the P107 audit without the P107.1 microplate protection policy")
    p107.add_argument("--disable-boundary-province-response", action="store_true",
                      help="run without the P107.3 terrain response to boundary provinces")
    p107.add_argument("--disable-p108-boundary-width-guard", action="store_true",
                      help="run without P108 boundary width guards")
    p107.add_argument("--disable-p108-high-mountain-coherence", action="store_true",
                      help="run without P108 high-mountain belt coherence response")
    p107.set_defaults(func=cmd_p107_audit)

    p110b = sub.add_parser(
        "p110b-seed-sweep",
        help="summarize existing P107 outputs as a P110B planform seed sweep",
    )
    p110b.add_argument(
        "inputs",
        nargs="+",
        help="P107 output directories, p107_audit_summary.json files, or terminal metrics",
    )
    p110b.add_argument("--out", default="out_p110b_seed_sweep")
    p110b.add_argument("--min-sample-size", type=int, default=5)
    p110b.add_argument("--max-soft-warning-rate", type=float, default=0.25)
    p110b.add_argument("--max-median-largest-share", type=float, default=0.58)
    p110b.add_argument("--max-p90-largest-share", type=float, default=0.62)
    p110b.set_defaults(func=cmd_p110b_seed_sweep)

    p107_render = sub.add_parser(
        "p107-render-arrays",
        help="render PNG visual QA assets from an existing P107 terminal array archive",
    )
    p107_render.add_argument(
        "input",
        help="P107 run directory or p107_terminal_metrics.json file",
    )
    p107_render.add_argument("--out", default="out_p107_array_render")
    p107_render.add_argument("--width", type=int, default=720)
    p107_render.add_argument("--height", type=int, default=360)
    p107_render.set_defaults(func=cmd_p107_render_arrays)

    ssr = sub.add_parser(
        "selected-snapshot-refine",
        help="derive a high-resolution selected-snapshot terrain refinement from P107 arrays",
    )
    ssr.add_argument(
        "input",
        help="P107 run directory or p107_terminal_metrics.json file",
    )
    ssr.add_argument("--out", default="out_selected_snapshot_refinement")
    ssr.add_argument("--target-cells", type=int, default=72000)
    ssr.add_argument("--width", type=int, default=1600)
    ssr.add_argument("--height", type=int, default=800)
    ssr.add_argument("--interpolation-k", type=int, default=6)
    ssr.add_argument("--detail-seed", type=int, default=72000)
    ssr.add_argument("--detail-strength", type=float, default=1.0)
    ssr.add_argument("--allow-process-islands", action="store_true")
    ssr.add_argument(
        "--render-groups",
        nargs="+",
        default=("all",),
        help=(
            "selected-snapshot QA render groups: all, p107, base, hydrology, "
            "marine, shelf, deep-ocean, submarine, island-atoll, coastal"
        ),
    )
    ssr.set_defaults(func=cmd_selected_snapshot_refine)

    ssr_render = sub.add_parser(
        "selected-snapshot-render-groups",
        help="render selected-snapshot QA groups from an existing refined output",
    )
    ssr_render.add_argument(
        "input",
        help="selected-snapshot output directory or selected_snapshot_refinement_metrics.json",
    )
    ssr_render.add_argument("--out", default=None)
    ssr_render.add_argument("--width", type=int, default=None)
    ssr_render.add_argument("--height", type=int, default=None)
    ssr_render.add_argument(
        "--render-groups",
        nargs="+",
        default=("all",),
        help=(
            "selected-snapshot QA render groups: all, p107, base, hydrology, "
            "marine, shelf, deep-ocean, submarine, island-atoll, coastal"
        ),
    )
    ssr_render.set_defaults(func=cmd_selected_snapshot_render_groups)

    p107_compare = sub.add_parser(
        "p107-compare",
        help="compare two existing P107 outputs for result-preserving optimization gates",
    )
    p107_compare.add_argument(
        "baseline",
        help="baseline P107 output root, run directory, or p107_terminal_metrics.json",
    )
    p107_compare.add_argument(
        "candidate",
        help="candidate P107 output root, run directory, or p107_terminal_metrics.json",
    )
    p107_compare.add_argument("--out", default="p107_equivalence_report.json")
    p107_compare.add_argument(
        "--metric-skip-key",
        action="append",
        default=[],
        help=(
            "top-level terminal metric key to ignore; defaults to profile and "
            "asset/archive metadata keys"
        ),
    )
    p107_compare.add_argument(
        "--float-atol",
        type=float,
        default=1e-12,
        help="absolute tolerance for numeric metric comparison",
    )
    p107_compare.set_defaults(func=cmd_p107_compare)

    tcb = sub.add_parser(
        "terminal-climate-biome",
        help="post-process accepted terminal terrains into climate and biome assets",
    )
    tcb.add_argument("--cells", type=int, default=8000)
    tcb.add_argument("--t-end", type=float, default=4500.0, dest="t_end")
    tcb.add_argument("--frames", type=int, default=4)
    tcb.add_argument("--max-workers", type=int, default=1)
    tcb.add_argument("--out", default="out_terminal_climate_biomes")
    tcb.add_argument("--no-render", action="store_true")
    tcb.add_argument(
        "--job",
        action="append",
        type=_parse_terminal_climate_job,
        default=[],
        metavar="PRESET:LABEL:SEED",
        help="override default six worlds; may be repeated",
    )
    tcb.set_defaults(func=cmd_terminal_climate_biome)

    ecr = sub.add_parser(
        "earth-climate-reference",
        help="build real-Earth reference manifest and same-grid calibration assets",
    )
    ecr.add_argument("--cells", nargs="+", type=int, default=[8000])
    ecr.add_argument("--width", type=int, default=720)
    ecr.add_argument("--height", type=int, default=360)
    ecr.add_argument("--out", default="out_earth_climate_reference")
    ecr.add_argument("--no-render", action="store_true")
    ecr.add_argument("--no-worldclim", action="store_true")
    ecr.add_argument("--no-koppen", action="store_true")
    ecr.add_argument("--no-noaa-psl", action="store_true")
    ecr.add_argument("--no-sst", action="store_true")
    ecr.add_argument("--no-ocean-currents", action="store_true")
    ecr.add_argument("--no-etopo2022", action="store_true")
    ecr.add_argument("--no-land-cover", action="store_true")
    ecr.add_argument(
        "--download-oscar",
        action="store_true",
        help="build the OSCAR 2001-2020 monthly current climatology cache if missing",
    )
    ecr.add_argument(
        "--download-land-cover",
        action="store_true",
        help="build the ESA CCI/C3S 2020 rendered-preview land-cover cache if missing",
    )
    ecr.add_argument("--climatology-start-year", type=int, default=1991)
    ecr.add_argument("--climatology-end-year", type=int, default=2020)
    ecr.add_argument(
        "--manifest-out",
        default=None,
        help="optional path to update the persistent source manifest",
    )
    ecr.set_defaults(func=cmd_earth_climate_reference)

    ecc = sub.add_parser(
        "earth-climate-compare",
        help="compare generated terminal climates against real-Earth reference metrics",
    )
    ecc.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecc.add_argument(
        "--terminal-summary",
        default="out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json",
    )
    ecc.add_argument("--out", default="out_earth_climate_comparison")
    ecc.add_argument("--no-contact-sheet", action="store_true")
    ecc.set_defaults(func=cmd_earth_climate_compare)

    ecf = sub.add_parser(
        "earth-climate-fit-report",
        help="summarize Earth comparison deltas into climate fitting priorities",
    )
    ecf.add_argument(
        "--comparison-summary",
        default=(
            "out_earth_climate_comparison_r4_20260705/"
            "earth_climate_comparison_summary.json"
        ),
    )
    ecf.add_argument("--out", default="out_earth_climate_fitting")
    ecf.add_argument(
        "--fail-on-guardrail",
        action="store_true",
        help="exit with status 2 when cross-preset guardrails fail",
    )
    ecf.set_defaults(func=cmd_earth_climate_fit_report)

    ecpg = sub.add_parser(
        "earth-climate-pattern-gate",
        help="evaluate Earthlike spatial climate-pattern envelopes",
    )
    ecpg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecpg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecpg.add_argument("--out", default="out_earth_climate_pattern_gate")
    ecpg.add_argument(
        "--fail-on-pattern",
        action="store_true",
        help="exit with status 2 when Earthlike pattern checks fail",
    )
    ecpg.set_defaults(func=cmd_earth_climate_pattern_gate)

    ecbg = sub.add_parser(
        "earth-climate-biome-gate",
        help="evaluate generated coarse biomes against Koppen/RESOLVE envelopes",
    )
    ecbg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecbg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecbg.add_argument("--out", default="out_earth_climate_biome_gate")
    ecbg.add_argument(
        "--fail-on-biome",
        action="store_true",
        help="exit with status 2 when Earthlike biome checks fail",
    )
    ecbg.set_defaults(func=cmd_earth_climate_biome_gate)

    ecsbg = sub.add_parser(
        "earth-climate-spatial-biome-gate",
        help="evaluate generated biome latitude organization against Earth envelopes",
    )
    ecsbg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecsbg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecsbg.add_argument("--out", default="out_earth_climate_spatial_biome_gate")
    ecsbg.add_argument(
        "--fail-on-spatial-biome",
        action="store_true",
        help="exit with status 2 when Earthlike spatial biome checks fail",
    )
    ecsbg.set_defaults(func=cmd_earth_climate_spatial_biome_gate)

    ecssg = sub.add_parser(
        "earth-climate-seasonal-subtype-gate",
        help="evaluate generated dry/wet-season subtype organization against Earth",
    )
    ecssg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecssg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecssg.add_argument("--out", default="out_earth_climate_seasonal_subtype_gate")
    ecssg.add_argument(
        "--fail-on-seasonal-subtype",
        action="store_true",
        help="exit with status 2 when Earthlike seasonal subtype checks fail",
    )
    ecssg.set_defaults(func=cmd_earth_climate_seasonal_subtype_gate)

    ecmzg = sub.add_parser(
        "earth-climate-mountain-zonation-gate",
        help="evaluate generated mountain ecological zonation against Earth envelopes",
    )
    ecmzg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecmzg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecmzg.add_argument("--out", default="out_earth_climate_mountain_zonation_gate")
    ecmzg.add_argument(
        "--fail-on-mountain-zonation",
        action="store_true",
        help="exit with status 2 when Earthlike mountain zonation checks fail",
    )
    ecmzg.set_defaults(func=cmd_earth_climate_mountain_zonation_gate)

    ecwlg = sub.add_parser(
        "earth-climate-windward-leeward-gate",
        help="evaluate generated windward/leeward mountain precipitation contrast against Earth",
    )
    ecwlg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecwlg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecwlg.add_argument("--out", default="out_earth_climate_windward_leeward_gate")
    ecwlg.add_argument(
        "--fail-on-windward-leeward",
        action="store_true",
        help="exit with status 2 when Earthlike windward/leeward checks fail",
    )
    ecwlg.set_defaults(func=cmd_earth_climate_windward_leeward_gate)

    ecmmsg = sub.add_parser(
        "earth-climate-monsoon-moisture-gate",
        help="evaluate C4a seasonal pressure, moisture access, and monsoon potential",
    )
    ecmmsg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r4_20260705/earth_reference_8000cells.npz",
    )
    ecmmsg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_f5wind3_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecmmsg.add_argument("--out", default="out_earth_climate_monsoon_moisture_gate")
    ecmmsg.add_argument(
        "--fail-on-monsoon-moisture",
        action="store_true",
        help="exit with status 2 when C4a monsoon/moisture checks fail",
    )
    ecmmsg.set_defaults(func=cmd_earth_climate_monsoon_moisture_gate)

    ecclg = sub.add_parser(
        "earth-climate-circulation-layout-gate",
        help="evaluate F2/F3 wind and ocean-current geography coupling",
    )
    ecclg.add_argument(
        "--earth-reference",
        default=(
            "out_earth_climate_reference_r6_landcover_20260706/"
            "earth_reference_8000cells.npz"
        ),
    )
    ecclg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c4k2_receiver_catchment_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecclg.add_argument("--out", default="out_earth_climate_circulation_layout_gate")
    ecclg.add_argument(
        "--fail-on-circulation-layout",
        action="store_true",
        help="exit with status 2 when F2/F3 circulation layout checks fail",
    )
    ecclg.set_defaults(func=cmd_earth_climate_circulation_layout_gate)

    rewr = sub.add_parser(
        "real-earth-wind-replay",
        help="compare one Earth-only R2 seasonal wind replay subgraph against reference winds",
    )
    rewr.add_argument(
        "--earth-reference",
        default=(
            "out_earth_climate_reference_r6_landcover_20260706/"
            "earth_reference_8000cells.npz"
        ),
    )
    rewr.add_argument(
        "--replay-arrays",
        default=(
            "out_real_earth_climate_replay_replay_r2_r4_phase4_render_20260706/"
            "terminal_climate_arrays.npz"
        ),
    )
    rewr.add_argument("--out", default="out_real_earth_wind_replay")
    rewr.set_defaults(func=cmd_real_earth_wind_replay)

    ecosg = sub.add_parser(
        "earth-climate-ocean-spatial-gate",
        help="evaluate F2 ocean-current and SST spatial structure",
    )
    ecosg.add_argument(
        "--earth-reference",
        default=(
            "out_earth_climate_reference_r6_landcover_20260706/"
            "earth_reference_8000cells.npz"
        ),
    )
    ecosg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c5a2_render_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecosg.add_argument("--out", default="out_earth_climate_ocean_spatial_gate")
    ecosg.add_argument(
        "--fail-on-ocean-spatial",
        action="store_true",
        help="exit with status 2 when F2 ocean spatial checks fail",
    )
    ecosg.set_defaults(func=cmd_earth_climate_ocean_spatial_gate)

    ecccg = sub.add_parser(
        "earth-climate-coupled-consistency-gate",
        help="evaluate F3 pressure/wind/moisture coupled consistency",
    )
    ecccg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r6_landcover_20260706/earth_reference_8000cells.npz",
    )
    ecccg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c5b1_ocean_spatial_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecccg.add_argument("--out", default="out_earth_climate_coupled_consistency_gate")
    ecccg.add_argument(
        "--fail-on-coupled-consistency",
        action="store_true",
        help="exit with status 2 when F3 coupled consistency checks fail",
    )
    ecccg.set_defaults(func=cmd_earth_climate_coupled_consistency_gate)

    ecshpg = sub.add_parser(
        "earth-climate-seasonal-hydro-placement-gate",
        help="evaluate F4 seasonal hydroclimate process placement",
    )
    ecshpg.add_argument(
        "--earth-reference",
        default="out_earth_climate_reference_r6_landcover_20260706/earth_reference_8000cells.npz",
    )
    ecshpg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c5c2_coupled_consistency_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecshpg.add_argument(
        "--out",
        default="out_earth_climate_seasonal_hydro_placement_gate",
    )
    ecshpg.add_argument(
        "--fail-on-seasonal-hydro-placement",
        action="store_true",
        help="exit with status 2 when F4 seasonal hydroclimate placement checks fail",
    )
    ecshpg.set_defaults(func=cmd_earth_climate_seasonal_hydro_placement_gate)

    eccvg = sub.add_parser(
        "earth-climate-coupling-convergence-gate",
        help="evaluate C5e bounded coupling feedback convergence",
    )
    eccvg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c5e1_hydro_feedback_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    eccvg.add_argument(
        "--out",
        default="out_earth_climate_coupling_convergence_gate",
    )
    eccvg.add_argument(
        "--fail-on-coupling-convergence",
        action="store_true",
        help="exit with status 2 when C5e coupling convergence checks fail",
    )
    eccvg.set_defaults(func=cmd_earth_climate_coupling_convergence_gate)

    echrg = sub.add_parser(
        "earth-climate-hydro-region-gate",
        help="evaluate C4d seasonal hydroclimate region objects",
    )
    echrg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c4d3_objects_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    echrg.add_argument("--out", default="out_earth_climate_hydro_region_gate")
    echrg.add_argument(
        "--fail-on-hydro-regions",
        action="store_true",
        help="exit with status 2 when C4d hydroclimate-region checks fail",
    )
    echrg.add_argument("--no-contact-sheet", action="store_true")
    echrg.set_defaults(func=cmd_earth_climate_hydro_region_gate)

    ecmfg = sub.add_parser(
        "earth-climate-moisture-flow-gate",
        help="evaluate C4e seasonal moisture-flow-network objects",
    )
    ecmfg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c4e1_flow_20260705/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecmfg.add_argument("--out", default="out_earth_climate_moisture_flow_gate")
    ecmfg.add_argument(
        "--fail-on-moisture-flow",
        action="store_true",
        help="exit with status 2 when C4e moisture-flow-network checks fail",
    )
    ecmfg.add_argument("--no-contact-sheet", action="store_true")
    ecmfg.set_defaults(func=cmd_earth_climate_moisture_flow_gate)

    ecmrg = sub.add_parser(
        "earth-climate-moisture-response-gate",
        help="evaluate C4f/C4j moisture-flow precipitation response, local budgets, and response objects",
    )
    ecmrg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c4j1_precip_objects_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecmrg.add_argument("--out", default="out_earth_climate_moisture_response_gate")
    ecmrg.add_argument(
        "--fail-on-moisture-response",
        action="store_true",
        help="exit with status 2 when C4f/C4j moisture-response checks fail",
    )
    ecmrg.add_argument("--no-contact-sheet", action="store_true")
    ecmrg.set_defaults(func=cmd_earth_climate_moisture_response_gate)

    ecrcg = sub.add_parser(
        "earth-climate-receiver-catchment-gate",
        help="evaluate C4k receiver catchments for source-basin/budget semantics",
    )
    ecrcg.add_argument(
        "--terminal-summary",
        default=(
            "out_terminal_climate_replay_c4j1_precip_objects_20260706/"
            "terminal_climate_replay_summary.json"
        ),
    )
    ecrcg.add_argument("--out", default="out_earth_climate_receiver_catchment_gate")
    ecrcg.add_argument(
        "--fail-on-receiver-catchment",
        action="store_true",
        help="exit with status 2 when C4k receiver-catchment checks fail",
    )
    ecrcg.set_defaults(func=cmd_earth_climate_receiver_catchment_gate)

    tcr = sub.add_parser(
        "terminal-climate-replay",
        help="rerun climate and static biomes on frozen terminal terrain arrays",
    )
    tcr.add_argument(
        "--terminal-summary",
        default="out_terminal_climate_biomes_20260705/terminal_climate_biome_summary.json",
    )
    tcr.add_argument(
        "--label",
        action="append",
        default=[],
        help="terminal world label to replay, e.g. earthlike_seed42; may be repeated",
    )
    tcr.add_argument("--out", default="out_terminal_climate_replay")
    tcr.add_argument("--no-render", action="store_true")
    tcr.set_defaults(func=cmd_terminal_climate_replay)

    recr = sub.add_parser(
        "real-earth-climate-replay",
        help="run current climate and biome layers on frozen real-Earth topography",
    )
    recr.add_argument(
        "--earth-reference",
        default=(
            "out_earth_climate_reference_r6_landcover_20260706/"
            "earth_reference_8000cells.npz"
        ),
    )
    recr.add_argument("--out", default="out_real_earth_climate_replay")
    recr.add_argument("--preset", default="earthlike", choices=sorted(PRESETS))
    recr.add_argument("--seed", type=int, default=20260706)
    recr.add_argument("--no-render", action="store_true")
    recr.set_defaults(func=cmd_real_earth_climate_replay)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
