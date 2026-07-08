import json

import numpy as np

from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS
from aevum.diagnostics.selected_snapshot_refinement import (
    _CoastalMorphology,
    _MarineMicrogeomorphology,
    _Microgeomorphology,
    SelectedSnapshotRefinementConfig,
    _block_local_neighborhood,
    _bridge_sparse_line_mask,
    _coastal_process_linework_zoom_centers,
    _fluvial_lacustrine_zoom_centers,
    _hydrology_zoom_centers,
    _marine_island_atoll_microrelief_delta,
    _marine_island_candidates,
    _marine_process_island_promotion,
    _marine_submarine_highland_morphology,
    _offshore_atoll_candidate_gate,
    _plateau_escarpment_segment_rank,
    _seamount_chain_peak_and_axis_rank,
    _shelf_slope_axis_rank,
    _spaced_rank_seed,
    _submarine_highlands_zoom_centers,
    _unsupported_inland_shallow_seaway_tuning,
    refine_selected_snapshot,
    render_selected_snapshot_refinement_assets,
)
from aevum.modules.terrain import (
    CONT_DETAIL_BASIN,
    CONT_DETAIL_OROGEN,
    CONT_DETAIL_PLATFORM,
    OCEAN_DEPTH_ABYSS,
    OCEAN_DEPTH_RIDGE,
    OCEAN_DEPTH_RESTRICTED,
)


def test_block_local_neighborhood_only_blocks_nearby_cells():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    center = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    blocked = np.zeros(grid.n, dtype=bool)

    _block_local_neighborhood(grid, blocked, center, passes=1)

    assert blocked[center]
    assert np.all(blocked[grid.neighbors[center]])
    assert np.count_nonzero(blocked) < 16


def test_bridge_sparse_line_mask_fills_single_cell_gap():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    start = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    middle = int(grid.neighbors[start][0])
    end_candidates = [int(c) for c in grid.neighbors[middle] if int(c) != start]
    end = end_candidates[0]
    seed = np.zeros(grid.n, dtype=bool)
    seed[start] = True
    seed[end] = True

    bridged = _bridge_sparse_line_mask(
        grid,
        seed,
        domain=np.ones(grid.n, dtype=bool),
        passes=1,
    )

    assert bridged[start]
    assert bridged[middle]
    assert bridged[end]
    assert np.count_nonzero(bridged) < 8


def test_shelf_slope_axis_rank_bridges_supported_single_cell_gap():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    start = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    middle = int(grid.neighbors[start][0])
    end_candidates = [int(c) for c in grid.neighbors[middle] if int(c) != start]
    end = end_candidates[0]
    weak = int(grid.neighbors[start][1])

    domain = np.zeros(grid.n, dtype=bool)
    domain[[start, middle, end, weak]] = True
    score = np.zeros(grid.n, dtype=np.float64)
    shelf_break = np.zeros(grid.n, dtype=np.float64)
    score[[start, end]] = 0.92
    score[middle] = 0.66
    score[weak] = 0.18
    shelf_break[domain] = 0.88

    rank = _shelf_slope_axis_rank(grid, domain, score, shelf_break)

    assert rank[start] > 0.0
    assert rank[middle] > 0.0
    assert rank[end] > 0.0
    assert rank[weak] == 0.0


def test_spaced_rank_seed_skips_immediate_neighbor_candidate():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    center = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    neighbor = int(grid.neighbors[center][0])
    far = int(np.argmin((grid.lat - 28.0) ** 2 + (grid.lon - 92.0) ** 2))
    values = np.zeros(grid.n, dtype=np.float64)
    values[center] = 1.0
    values[neighbor] = 0.98
    values[far] = 0.82

    seed = _spaced_rank_seed(
        grid,
        np.asarray([center, neighbor, far], dtype=np.int64),
        values,
        max_count=3,
        spacing_passes=1,
    )

    assert seed[center]
    assert not seed[neighbor]
    assert seed[far]


def test_seamount_chain_peak_axis_keeps_peaks_spaced_and_bridges_support():
    grid = SphereGrid.fibonacci(320, CONSTANTS.EARTH_RADIUS)
    start = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    middle = int(grid.neighbors[start][0])
    end_candidates = [int(c) for c in grid.neighbors[middle] if int(c) != start]
    end = end_candidates[0]
    near_start = int(grid.neighbors[start][1])
    filler = [
        int(c)
        for c in np.flatnonzero(np.ones(grid.n, dtype=bool))
        if int(c) not in {start, middle, end, near_start}
    ][:6]

    domain = np.zeros(grid.n, dtype=bool)
    domain[[start, middle, end, near_start]] = True
    domain[filler] = True
    score = np.zeros(grid.n, dtype=np.float64)
    score[start] = 0.95
    score[near_start] = 0.92
    score[middle] = 0.66
    score[end] = 0.90
    score[filler] = 0.50

    peak_rank, axis_rank = _seamount_chain_peak_and_axis_rank(grid, domain, score)

    assert peak_rank[start] > 0.0
    assert peak_rank[near_start] == 0.0
    assert peak_rank[end] > 0.0
    assert axis_rank[middle] > 0.0


def test_plateau_escarpment_segment_rank_does_not_force_whole_edge():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    plateau = (
        (np.abs(grid.lat) < 18.0)
        & (grid.lon > -42.0)
        & (grid.lon < 42.0)
    )
    edge = np.zeros(grid.n, dtype=bool)
    for cell in np.flatnonzero(plateau):
        nbs = grid.neighbors[int(cell)]
        if nbs.size and np.any(~plateau[nbs]):
            edge[int(cell)] = True
    edge_cells = np.flatnonzero(edge)
    assert edge_cells.size > 8
    ridge_texture = np.zeros(grid.n, dtype=np.float64)
    ridge_texture[edge_cells[: max(3, edge_cells.size // 4)]] = 1.0

    rank = _plateau_escarpment_segment_rank(
        grid,
        edge,
        np.ones(grid.n, dtype=bool),
        plateau,
        np.zeros(grid.n, dtype=bool),
        ridge_texture,
        np.zeros(grid.n, dtype=np.float64),
    )

    active_edge = edge & (rank > 0.0)
    assert np.count_nonzero(active_edge) > 0
    assert np.count_nonzero(active_edge) < edge_cells.size
    assert np.max(rank[edge]) >= 0.82


def test_submarine_highland_microcontinent_edge_rank_does_not_fill_interior():
    grid = SphereGrid.fibonacci(1800, CONSTANTS.EARTH_RADIUS)
    microcontinent = (
        (np.abs(grid.lat) < 17.0)
        & (grid.lon > -36.0)
        & (grid.lon < 36.0)
    )
    edge = np.zeros(grid.n, dtype=bool)
    for cell in np.flatnonzero(microcontinent):
        nbs = grid.neighbors[int(cell)]
        if nbs.size and np.any(~microcontinent[nbs]):
            edge[int(cell)] = True
    edge_halo = edge.copy()
    for cell in np.flatnonzero(edge):
        edge_halo[grid.neighbors[int(cell)]] = True
    interior = microcontinent & ~edge_halo
    assert np.count_nonzero(edge) > 0
    assert np.count_nonzero(interior) > 0

    rel = np.full(grid.n, -1800.0, dtype=np.float64)
    potential_rel = rel + microcontinent.astype(np.float64) * 420.0
    zeros = np.zeros(grid.n, dtype=bool)
    delta, _, _, plateau_edge_rank, _ = _marine_submarine_highland_morphology(
        grid,
        rel,
        potential_rel,
        zeros,
        zeros,
        microcontinent,
        zeros,
        zeros,
        zeros,
        np.full(grid.n, OCEAN_DEPTH_ABYSS, dtype=np.int32),
        np.full(grid.n, 8.0, dtype=np.float64),
        np.full(grid.n, 0.5, dtype=np.float64),
        np.zeros(grid.n, dtype=np.float64),
        np.ones(grid.n, dtype=bool),
    )

    assert np.max(plateau_edge_rank[edge]) >= 0.72
    assert np.max(plateau_edge_rank[interior]) == 0.0
    assert np.max(np.abs(delta[interior])) == 0.0


def test_submarine_highland_zoom_prefers_oceanic_seamount_peak():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    near = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    far = int(np.argmin((grid.lat - 20.0) ** 2 + (grid.lon - 100.0) ** 2))
    zeros = np.zeros(grid.n, dtype=np.float64)
    false = np.zeros(grid.n, dtype=bool)
    seamount_peak = zeros.copy()
    seamount_peak[near] = 1.0
    seamount_peak[far] = 0.82
    coast_distance = np.full(grid.n, 8.0, dtype=np.float64)
    coast_distance[near] = 1.0
    coast_distance[far] = 12.0
    refined = np.full(grid.n, -1800.0, dtype=np.float64)
    refined[near] = -90.0
    refined[far] = -2300.0
    marine = _MarineMicrogeomorphology(
        refined_rel=refined,
        marine_delta=zeros,
        ocean_coast_distance=coast_distance,
        shelf_break_rank=zeros,
        shelf_slope_microrelief_delta=zeros,
        deep_ocean_fabric_delta=zeros,
        fracture_zone_rank=zeros,
        abyssal_plain_fabric_rank=zeros,
        island_candidate_rank=zeros,
        atoll_candidate_rank=zeros,
        reef_atoll_mask=false,
        marine_shoal_mask=false,
        seamount_shoal_mask=false,
        oceanic_plateau_shoal_mask=false,
        microcontinent_shoal_mask=false,
        island_arc_shoal_mask=false,
        process_island_candidate_mask=false,
        atoll_candidate_mask=false,
        reef_rim_rank=zeros,
        atoll_lagoon_rank=zeros,
        fringing_reef_rank=zeros,
        process_island_promotion_rank=zeros,
        process_island_promotion_delta=zeros,
        process_island_promoted_mask=false,
        atoll_islet_promoted_mask=false,
        islet_microshape_rank=zeros,
        atoll_microshape_rank=zeros,
        island_atoll_microrelief_delta=zeros,
        submarine_highland_delta=zeros,
        inland_seaway_tuning_delta=zeros,
        inland_seaway_tuning_mask=false,
        inland_seaway_tuning_rank=zeros,
        inland_seaway_landback_mask=false,
        seamount_peak_rank=seamount_peak,
        seamount_apron_rank=zeros,
        oceanic_plateau_edge_rank=zeros,
        abyssal_hill_field_rank=zeros,
    )

    label, lon, lat = _submarine_highlands_zoom_centers(grid, marine)[0]

    assert label == "oceanic seamount peak"
    assert abs(lon - float(grid.lon[far])) < 1.0e-9
    assert abs(lat - float(grid.lat[far])) < 1.0e-9


def test_coastal_process_zoom_centers_avoid_duplicate_local_windows():
    grid = SphereGrid.fibonacci(900, CONSTANTS.EARTH_RADIUS)
    shared = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    alternate_estuary = int(
        np.argmin((grid.lat - 24.0) ** 2 + (grid.lon - 112.0) ** 2)
    )
    barrier = int(np.argmin((grid.lat + 26.0) ** 2 + (grid.lon + 118.0) ** 2))
    zeros = np.zeros(grid.n, dtype=np.float64)
    delta = zeros.copy()
    estuary = zeros.copy()
    spit = zeros.copy()
    delta[shared] = 1.0
    estuary[shared] = 0.98
    estuary[alternate_estuary] = 0.82
    spit[barrier] = 0.75
    coastal = _CoastalMorphology(
        refined_rel=zeros,
        coastal_delta=zeros,
        coastal_process_microrelief_delta=zeros,
        coastal_depositional_microrelief_delta=zeros,
        coastal_plain_rank=zeros,
        coastal_cliff_rank=zeros,
        shoreface_rank=zeros,
        barrier_lagoon_rank=zeros,
        estuary_rank=zeros,
        delta_distributary_rank=delta,
        estuary_funnel_rank=estuary,
        barrier_spit_rank=spit,
        delta_mouth_bar_rank=zeros,
        estuary_tidal_channel_rank=zeros,
        coastal_depositional_plain_rank=zeros,
        strandplain_rank=zeros,
        tidal_flat_rank=zeros,
    )

    centers = _coastal_process_linework_zoom_centers(grid, coastal)

    assert centers[0][0] == "delta distributaries"
    assert centers[1][0] == "estuary funnel"
    assert abs(centers[0][1] - float(grid.lon[shared])) < 1.0e-9
    assert abs(centers[0][2] - float(grid.lat[shared])) < 1.0e-9
    assert abs(centers[1][1] - float(grid.lon[alternate_estuary])) < 1.0e-9
    assert abs(centers[1][2] - float(grid.lat[alternate_estuary])) < 1.0e-9


def _blank_microgeomorphology(grid: SphereGrid) -> _Microgeomorphology:
    zeros = np.zeros(grid.n, dtype=np.float64)
    false = np.zeros(grid.n, dtype=bool)
    return _Microgeomorphology(
        refined_rel=zeros.copy(),
        hydrology_delta=zeros.copy(),
        fluvial_microrelief_delta=zeros.copy(),
        lowland_alluvial_microrelief_delta=zeros.copy(),
        flow_accumulation=zeros.copy(),
        river_rank=zeros.copy(),
        river_path_rank=zeros.copy(),
        basin_trunk_rank=zeros.copy(),
        floodplain_rank=zeros.copy(),
        meander_belt_rank=zeros.copy(),
        meander_scroll_rank=zeros.copy(),
        floodplain_swale_rank=zeros.copy(),
        alluvial_fan_rank=zeros.copy(),
        lowland_plain_rank=zeros.copy(),
        piedmont_apron_rank=zeros.copy(),
        lake_basin_rank=zeros.copy(),
        lake_shoreline_rank=zeros.copy(),
        delta_fan_rank=zeros.copy(),
        drainage_basin_id=np.zeros(grid.n, dtype=np.int32),
        river_receiver=np.full(grid.n, -1, dtype=np.int64),
        lake_mask=false.copy(),
        delta_mask=false.copy(),
        delta_plain_mask=false.copy(),
        land_coast_distance=np.full(grid.n, -1.0, dtype=np.float64),
    )


def test_hydrology_zoom_centers_choose_landward_river_mouth():
    grid = SphereGrid.fibonacci(540, CONSTANTS.EARTH_RADIUS)
    mouth = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    offshore = int(grid.neighbors[mouth][0])
    micro = _blank_microgeomorphology(grid)
    rel = np.full(grid.n, 180.0, dtype=np.float64)
    rel[offshore] = -90.0
    micro.land_coast_distance[mouth] = 0.0
    micro.river_path_rank[mouth] = 0.78
    micro.river_rank[mouth] = 0.64
    micro.delta_plain_mask[mouth] = True
    micro.delta_fan_rank[offshore] = 1.0

    centers = _hydrology_zoom_centers(grid, rel, micro)
    mouth_center = next(center for center in centers if center[0] == "river mouth")

    assert abs(mouth_center[1] - float(grid.lon[mouth])) < 1.0e-9
    assert abs(mouth_center[2] - float(grid.lat[mouth])) < 1.0e-9


def test_fluvial_zoom_centers_prefer_temperate_lake_shoreline_candidate():
    grid = SphereGrid.fibonacci(900, CONSTANTS.EARTH_RADIUS)
    polar_lake = int(np.argmin((grid.lat - 68.0) ** 2 + (grid.lon - 30.0) ** 2))
    temperate_lake = int(np.argmin((grid.lat - 32.0) ** 2 + (grid.lon + 84.0) ** 2))
    micro = _blank_microgeomorphology(grid)
    micro.lake_shoreline_rank[polar_lake] = 1.0
    micro.lake_shoreline_rank[temperate_lake] = 0.82

    centers = _fluvial_lacustrine_zoom_centers(grid, micro)
    lake_center = next(center for center in centers if center[0] == "lake shoreline")

    assert abs(lake_center[1] - float(grid.lon[temperate_lake])) < 1.0e-9
    assert abs(lake_center[2] - float(grid.lat[temperate_lake])) < 1.0e-9


def test_offshore_atoll_gate_rejects_nearshore_arc_and_accepts_remote_highland():
    ocean_distance = np.asarray([1.0, 3.0, 4.0, 5.0], dtype=np.float64)
    false = np.zeros(4, dtype=bool)
    seamount = false.copy()
    arc = false.copy()
    microcontinent = false.copy()
    seamount[1] = True
    arc[0] = True
    arc[2] = True
    microcontinent[3] = True

    gate = _offshore_atoll_candidate_gate(
        ocean_distance,
        seamount=seamount,
        plateau=false,
        microcontinent=microcontinent,
        island_arc=arc,
    )

    assert not gate[0]
    assert gate[1]
    assert gate[2]
    assert gate[3]


def test_marine_island_candidates_penalize_nearshore_microcontinent():
    grid = SphereGrid.fibonacci(720, CONSTANTS.EARTH_RADIUS)
    near = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    far = int(np.argmin((grid.lat - 12.0) ** 2 + (grid.lon - 96.0) ** 2))
    rel = np.full(grid.n, -760.0, dtype=np.float64)
    potential = np.full(grid.n, -180.0, dtype=np.float64)
    process_high = np.zeros(grid.n, dtype=bool)
    microcontinent = np.zeros(grid.n, dtype=bool)
    process_high[[near, far]] = True
    microcontinent[[near, far]] = True
    ocean_distance = np.full(grid.n, 6.0, dtype=np.float64)
    ocean_distance[near] = 1.0
    ridge_texture = np.full(grid.n, 0.65, dtype=np.float64)

    rank, mask = _marine_island_candidates(
        grid,
        rel,
        potential,
        process_high,
        np.zeros(grid.n, dtype=bool),
        np.zeros(grid.n, dtype=bool),
        np.zeros(grid.n, dtype=bool),
        microcontinent,
        np.zeros(grid.n, dtype=bool),
        ridge_texture,
        np.ones(grid.n, dtype=bool),
        ocean_distance,
    )

    assert rank[near] == 0.0
    assert not mask[near]
    assert rank[far] > 0.0
    assert mask[far]


def test_process_island_promotion_classifies_hybrid_rim_candidate_as_atoll_islet():
    grid = SphereGrid.fibonacci(180, CONSTANTS.EARTH_RADIUS)
    cell = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    zeros = np.zeros(grid.n, dtype=np.float64)
    false = np.zeros(grid.n, dtype=bool)
    rel = np.full(grid.n, -90.0, dtype=np.float64)
    potential = np.full(grid.n, -45.0, dtype=np.float64)
    island_rank = zeros.copy()
    atoll_rank = zeros.copy()
    rim_rank = zeros.copy()
    island_rank[cell] = 0.85
    atoll_rank[cell] = 0.80
    rim_rank[cell] = 0.70
    process_candidate = false.copy()
    atoll_candidate = false.copy()
    microcontinent = false.copy()
    process_candidate[cell] = True
    atoll_candidate[cell] = True
    microcontinent[cell] = True

    _, promotion_rank, promoted, atoll_promoted = _marine_process_island_promotion(
        grid,
        rel,
        potential,
        island_rank,
        atoll_rank,
        rim_rank,
        zeros.copy(),
        zeros.copy(),
        process_candidate,
        atoll_candidate,
        false.copy(),
        false.copy(),
        microcontinent,
        false.copy(),
        false.copy(),
        np.ones(grid.n, dtype=bool),
    )

    assert promoted[cell]
    assert atoll_promoted[cell]
    assert promotion_rank[cell] > 0.0


def test_unsupported_inland_shallow_seaway_tuning_narrows_gateway_only():
    grid = SphereGrid.fibonacci(1200, CONSTANTS.EARTH_RADIUS)
    unsupported = (
        (np.abs(grid.lat) < 4.5)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )
    supported = (
        (np.abs(grid.lat - 22.0) < 4.5)
        & (grid.lon > -95.0)
        & (grid.lon < 95.0)
    )
    ocean = unsupported | supported
    rel = np.where(ocean, -220.0, 420.0).astype(np.float64)
    distance = np.where(ocean, 1.4, 99.0).astype(np.float64)
    depth = np.zeros(grid.n, dtype=np.int32)
    depth[ocean] = OCEAN_DEPTH_RESTRICTED
    shelf_width = np.where(ocean, 1.0, 0.0).astype(np.float64)
    gateway_id = np.full(grid.n, -1, dtype=np.int32)
    gateway_id[ocean] = 0
    gateway_system_id = gateway_id.copy()
    process_support = supported.copy()

    delta, tuned, rank, landback = _unsupported_inland_shallow_seaway_tuning(
        grid,
        rel,
        ocean,
        distance,
        depth,
        shelf_width,
        gateway_id,
        gateway_system_id,
        process_support,
    )

    assert np.count_nonzero(tuned & unsupported) > 0
    assert np.count_nonzero(landback & unsupported) > 0
    assert np.any(delta[tuned & ~landback] < 0.0)
    assert np.all(delta[landback] > 0.0)
    assert np.max(rank[tuned & unsupported]) > 0.0
    assert np.count_nonzero(tuned & supported) == 0
    assert np.count_nonzero(landback & supported) == 0


def test_island_atoll_microrelief_keeps_halo_weaker_than_core():
    grid = SphereGrid.fibonacci(420, CONSTANTS.EARTH_RADIUS)
    center = int(np.argmin(grid.lat * grid.lat + grid.lon * grid.lon))
    halo = int(grid.neighbors[center][0])

    rel = np.full(grid.n, -180.0, dtype=np.float64)
    refined = np.full(grid.n, -92.0, dtype=np.float64)
    islet_rank = np.zeros(grid.n, dtype=np.float64)
    islet_rank[center] = 1.0
    reef_rank = np.zeros(grid.n, dtype=np.float64)
    reef_rank[center] = 0.7
    reef_rank[halo] = 0.7

    delta = _marine_island_atoll_microrelief_delta(
        grid,
        rel,
        refined,
        islet_rank,
        np.zeros(grid.n, dtype=np.float64),
        reef_rank,
        np.zeros(grid.n, dtype=np.float64),
        np.zeros(grid.n, dtype=np.float64),
        np.zeros(grid.n, dtype=bool),
        np.zeros(grid.n, dtype=bool),
        np.ones(grid.n, dtype=bool),
        np.full(grid.n, 0.5, dtype=np.float64),
        np.ones(grid.n, dtype=bool),
        1234,
    )

    assert delta[center] > 0.0
    assert delta[halo] < 0.45 * delta[center]


def test_selected_snapshot_refinement_writes_parented_qa_pack(tmp_path):
    grid = SphereGrid.fibonacci(180, CONSTANTS.EARTH_RADIUS)
    land = grid.lat > -8.0
    rel = np.where(land, 260.0 + 7.0 * grid.lat, -3600.0)
    rel[(grid.lat > 35.0) & (grid.lon > -30.0) & (grid.lon < 70.0)] = 1200.0
    rel[(grid.lat < -35.0) & ~land] = -5200.0
    detail = np.where(land, CONT_DETAIL_PLATFORM, 0).astype(np.int32)
    detail[(grid.lat > 35.0) & land] = CONT_DETAIL_OROGEN
    detail[(grid.lat > 0.0) & (grid.lat < 15.0) & land] = CONT_DETAIL_BASIN
    hierarchy = np.zeros(grid.n, dtype=np.int32)
    hierarchy[(grid.lat > 38.0) & land] = 2
    spine = np.zeros(grid.n, dtype=np.int32)
    spine[(grid.lat > 44.0) & (grid.lat < 50.0) & land] = 3
    ridge = ((np.abs(grid.lat + 25.0) < 6.0) & ~land).astype(np.uint8)
    trench = ((grid.lat < -42.0) & (grid.lon > 10.0)).astype(np.uint8)
    seamount = ((grid.lat < -15.0) & (grid.lat > -35.0) & ~land).astype(np.uint8)
    transform = (
        (np.abs(grid.lon + 45.0) < 7.0)
        & (grid.lat < -10.0)
        & (grid.lat > -55.0)
        & ~land
    ).astype(np.uint8)
    abyssal_plain = (
        (grid.lon < -35.0)
        & (grid.lat < -10.0)
        & (grid.lat > -45.0)
        & (np.abs(grid.lat + 25.0) > 9.0)
        & ~land
    ).astype(np.uint8)
    depth_province = np.where(land, 0, OCEAN_DEPTH_ABYSS).astype(np.int32)
    depth_province[ridge.astype(bool)] = OCEAN_DEPTH_RIDGE

    arrays_path = tmp_path / "p107_terminal_arrays.npz"
    np.savez_compressed(
        arrays_path,
        grid_lat=grid.lat.astype(np.float32),
        grid_lon=grid.lon.astype(np.float32),
        grid_cell_area_m2=grid.cell_area.astype(np.float64),
        field__terrain_elevation_m=rel.astype(np.float32),
        field__tectonics_plate_id=(np.arange(grid.n) % 8).astype(np.int32),
        field__crust_type=land.astype(np.int32),
        field__crust_age_myr=np.where(land, 1600.0, 155.0).astype(np.float32),
        field__ocean_depth_province=depth_province,
        field__terrain_continental_detail=detail,
        field__terrain_continental_detail_region_code=detail,
        field__terrain_orogenic_parent_hierarchy=hierarchy,
        field__terrain_orogenic_hierarchy_spine=spine,
        field__terrain_orogenic_shoulder_halo=(hierarchy > 0).astype(np.int32),
        field__terrain_orogenic_highland_apron=np.zeros(grid.n, dtype=np.int32),
        boundary__ridge=ridge,
        boundary__trench=trench,
        boundary__transform=transform,
        object__province_mid_ocean_ridge=ridge,
        object__margin_trench=trench,
        object__transform_fault=transform,
        object__seamount_chain=seamount,
        object__abyssal_plain=abyssal_plain,
        object__oceanic_plateau=np.zeros(grid.n, dtype=np.uint8),
        object__microcontinent=np.zeros(grid.n, dtype=np.uint8),
        object__island_arc=np.zeros(grid.n, dtype=np.uint8),
        object__mountain_orogen=(hierarchy > 0).astype(np.uint8),
        object__parent_orogen_crest=(spine > 0).astype(np.uint8),
    )
    metrics = {
        "schema": "aevum.p107_terminal_world_audit.v1",
        "cells": grid.n,
        "sea_level_m": 0.0,
        "array_archive": {
            "path": str(arrays_path),
            "manifest": {
                "fields": {
                    "terrain.elevation_m": "field__terrain_elevation_m",
                    "tectonics.plate_id": "field__tectonics_plate_id",
                    "crust.type": "field__crust_type",
                    "crust.age_myr": "field__crust_age_myr",
                    "ocean.depth_province": "field__ocean_depth_province",
                    "terrain.continental_detail": "field__terrain_continental_detail",
                    "terrain.continental_detail_region_code": (
                        "field__terrain_continental_detail_region_code"
                    ),
                    "terrain.orogenic_parent_hierarchy": (
                        "field__terrain_orogenic_parent_hierarchy"
                    ),
                    "terrain.orogenic_hierarchy_spine": (
                        "field__terrain_orogenic_hierarchy_spine"
                    ),
                    "terrain.orogenic_shoulder_halo": (
                        "field__terrain_orogenic_shoulder_halo"
                    ),
                    "terrain.orogenic_highland_apron": (
                        "field__terrain_orogenic_highland_apron"
                    ),
                }
            },
        },
        "p109_hypsometry_comparison": {"within_p109_envelope": True},
        "p110a_modern_planform": {
            "summary": {"land_fraction": float(np.mean(land))},
            "warning_flags": [],
        },
    }
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    summary = refine_selected_snapshot(
        SelectedSnapshotRefinementConfig(
            source=metrics_path,
            outdir=tmp_path / "refined",
            target_cells=540,
            width=96,
            height=48,
            detail_seed=1234,
        )
    )

    assert summary["schema"] == "aevum.selected_snapshot_refinement.v1"
    assert summary["source_cells"] == grid.n
    assert summary["target_cells"] == 540
    assert summary["land_ocean_sign_flip_fraction"] == 0.0
    assert summary["detail_delta_abs_p95_m"] > 0.0
    assert summary["marine_delta_nonzero_abs_p95_m"] > 0.0
    assert summary["river_cell_fraction_land"] > 0.0
    assert summary["drainage_basin_count"] >= 0
    assert summary["meander_belt_cell_fraction_land"] >= 0.0
    assert summary["meander_scroll_cell_fraction_land"] >= 0.0
    assert summary["floodplain_swale_cell_fraction_land"] >= 0.0
    assert summary["lake_basin_cell_fraction_land"] >= summary["lake_cell_fraction_land"]
    assert summary["lake_shoreline_cell_fraction_land"] >= summary["lake_cell_fraction_land"]
    assert summary["delta_fan_cell_fraction_ocean"] >= summary["delta_cell_fraction_ocean"]
    assert summary["process_island_candidate_cell_fraction_ocean"] >= 0.0
    assert summary["atoll_candidate_cell_fraction_ocean"] >= 0.0
    assert summary["reef_rim_cell_fraction_ocean"] >= 0.0
    assert summary["atoll_lagoon_cell_fraction_ocean"] >= 0.0
    assert summary["fringing_reef_cell_fraction_ocean"] >= 0.0
    assert summary["process_island_promotion_rank_cell_fraction_ocean"] == 0.0
    assert summary["process_island_promoted_cell_fraction_parent_ocean"] == 0.0
    assert summary["atoll_islet_promoted_cell_fraction_parent_ocean"] == 0.0
    assert summary["islet_microshape_rank_cell_fraction_ocean"] >= 0.0
    assert summary["atoll_microshape_rank_cell_fraction_ocean"] >= 0.0
    assert summary["island_atoll_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["island_atoll_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["island_atoll_microrelief_cell_fraction_ocean"] >= 0.0
    assert summary["seamount_peak_cell_fraction_ocean"] >= 0.0
    assert summary["seamount_apron_cell_fraction_ocean"] >= 0.0
    assert summary["oceanic_plateau_edge_cell_fraction_ocean"] >= 0.0
    assert summary["abyssal_hill_field_cell_fraction_ocean"] >= 0.0
    assert summary["inland_seaway_tuning_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["inland_seaway_tuning_cell_fraction_ocean"] >= 0.0
    assert summary["inland_seaway_landback_cell_fraction_parent_ocean"] >= 0.0
    assert summary["coastal_plain_cell_fraction_land"] >= 0.0
    assert summary["coastal_cliff_cell_fraction_land"] >= 0.0
    assert summary["shoreface_cell_fraction_ocean"] >= 0.0
    assert summary["barrier_lagoon_cell_fraction_ocean"] >= 0.0
    assert summary["estuary_cell_fraction_ocean"] >= 0.0
    assert summary["delta_distributary_cell_fraction_ocean"] >= 0.0
    assert summary["estuary_funnel_cell_fraction_ocean"] >= 0.0
    assert summary["barrier_spit_cell_fraction_ocean"] >= 0.0
    assert summary["delta_mouth_bar_cell_fraction_ocean"] >= 0.0
    assert summary["estuary_tidal_channel_cell_fraction_ocean"] >= 0.0
    assert summary["hydrology_delta_abs_p95_m"] > 0.0
    assert summary["fluvial_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["fluvial_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["fluvial_microrelief_cell_fraction_land"] >= 0.0
    assert summary["lowland_alluvial_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["lowland_alluvial_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["lowland_alluvial_microrelief_cell_fraction_land"] >= 0.0
    assert summary["alluvial_fan_rank_cell_fraction_land"] >= 0.0
    assert summary["lowland_plain_rank_cell_fraction_land"] >= 0.0
    assert summary["piedmont_apron_rank_cell_fraction_land"] >= 0.0
    assert summary["coastal_process_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["coastal_process_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["coastal_process_microrelief_cell_fraction_ocean"] >= 0.0
    assert summary["coastal_depositional_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["coastal_depositional_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["coastal_depositional_microrelief_cell_fraction_land"] >= 0.0
    assert summary["coastal_depositional_plain_rank_cell_fraction_land"] >= 0.0
    assert summary["strandplain_rank_cell_fraction_land"] >= 0.0
    assert summary["tidal_flat_rank_cell_fraction_land"] >= 0.0
    assert summary["shelf_slope_microrelief_delta_abs_p95_m"] >= 0.0
    assert summary["shelf_slope_microrelief_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["shelf_slope_microrelief_cell_fraction_ocean"] >= 0.0
    assert summary["deep_ocean_fabric_delta_abs_p95_m"] >= 0.0
    assert summary["deep_ocean_fabric_delta_nonzero_abs_p95_m"] >= 0.0
    assert summary["deep_ocean_fabric_cell_fraction_ocean"] >= 0.0
    assert summary["fracture_zone_rank_cell_fraction_ocean"] >= 0.0
    assert summary["abyssal_plain_fabric_rank_cell_fraction_ocean"] >= 0.0
    assert (tmp_path / "refined" / "selected_snapshot_refined_arrays.npz").exists()
    assert (tmp_path / "refined" / "rendered" / "elevation.png").exists()
    assert (
        tmp_path / "refined" / "rendered" / "parent_vs_refined_elevation.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" / "selected_snapshot_hydrology.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" / "selected_snapshot_hydrology_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" / "selected_snapshot_drainage_basins.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_fluvial_lacustrine_microshapes.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_fluvial_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_fluvial_microrelief_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_lowland_alluvial_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_marine_microgeomorphology.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" / "selected_snapshot_marine_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_marine_object_classes.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_shelf_slope_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_shelf_slope_microrelief_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_deep_ocean_fabric_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_deep_ocean_fabric_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_submarine_highlands.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_submarine_highlands_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_island_atoll_candidates.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_process_island_promotion.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_process_island_promotion_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_island_atoll_microshapes.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_island_atoll_microshapes_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_island_atoll_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_island_atoll_microrelief_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_reef_atoll_morphology.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_reef_atoll_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_morphology.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_process_linework.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_process_linework_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_process_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_process_microrelief_zoom_sheet.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_depositional_microrelief_delta.png"
    ).exists()
    assert (
        tmp_path / "refined" / "rendered" /
        "selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png"
    ).exists()
    with np.load(tmp_path / "refined" / "selected_snapshot_refined_arrays.npz") as z:
        assert "field__selected_snapshot_marine_delta_m" in z.files
        assert "field__selected_snapshot_ocean_coast_distance_passes" in z.files
        assert "field__selected_snapshot_shelf_break_rank" in z.files
        assert "field__selected_snapshot_shelf_slope_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_deep_ocean_fabric_delta_m" in z.files
        assert "field__selected_snapshot_fracture_zone_rank" in z.files
        assert "field__selected_snapshot_abyssal_plain_fabric_rank" in z.files
        assert "field__selected_snapshot_island_candidate_rank" in z.files
        assert "field__selected_snapshot_atoll_candidate_rank" in z.files
        assert "field__selected_snapshot_reef_rim_rank" in z.files
        assert "field__selected_snapshot_atoll_lagoon_rank" in z.files
        assert "field__selected_snapshot_fringing_reef_rank" in z.files
        assert "field__selected_snapshot_process_island_promotion_rank" in z.files
        assert "field__selected_snapshot_process_island_promotion_delta_m" in z.files
        assert "field__selected_snapshot_islet_microshape_rank" in z.files
        assert "field__selected_snapshot_atoll_microshape_rank" in z.files
        assert "field__selected_snapshot_island_atoll_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_submarine_highland_delta_m" in z.files
        assert "field__selected_snapshot_inland_seaway_tuning_delta_m" in z.files
        assert "field__selected_snapshot_inland_seaway_tuning_rank" in z.files
        assert "mask__selected_snapshot_inland_seaway_landback" in z.files
        assert "field__selected_snapshot_seamount_peak_rank" in z.files
        assert "field__selected_snapshot_seamount_apron_rank" in z.files
        assert "field__selected_snapshot_oceanic_plateau_edge_rank" in z.files
        assert "field__selected_snapshot_abyssal_hill_field_rank" in z.files
        assert "field__selected_snapshot_fluvial_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_lowland_alluvial_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_alluvial_fan_rank" in z.files
        assert "field__selected_snapshot_lowland_plain_rank" in z.files
        assert "field__selected_snapshot_piedmont_apron_rank" in z.files
        assert "field__selected_snapshot_coastal_delta_m" in z.files
        assert "field__selected_snapshot_coastal_process_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_coastal_depositional_microrelief_delta_m" in z.files
        assert "field__selected_snapshot_coastal_plain_rank" in z.files
        assert "field__selected_snapshot_coastal_cliff_rank" in z.files
        assert "field__selected_snapshot_shoreface_rank" in z.files
        assert "field__selected_snapshot_barrier_lagoon_rank" in z.files
        assert "field__selected_snapshot_estuary_rank" in z.files
        assert "field__selected_snapshot_delta_distributary_rank" in z.files
        assert "field__selected_snapshot_estuary_funnel_rank" in z.files
        assert "field__selected_snapshot_barrier_spit_rank" in z.files
        assert "field__selected_snapshot_delta_mouth_bar_rank" in z.files
        assert "field__selected_snapshot_estuary_tidal_channel_rank" in z.files
        assert "field__selected_snapshot_coastal_depositional_plain_rank" in z.files
        assert "field__selected_snapshot_strandplain_rank" in z.files
        assert "field__selected_snapshot_tidal_flat_rank" in z.files
        assert "field__selected_snapshot_basin_trunk_rank" in z.files
        assert "field__selected_snapshot_floodplain_rank" in z.files
        assert "field__selected_snapshot_meander_belt_rank" in z.files
        assert "field__selected_snapshot_meander_scroll_rank" in z.files
        assert "field__selected_snapshot_floodplain_swale_rank" in z.files
        assert "field__selected_snapshot_lake_basin_rank" in z.files
        assert "field__selected_snapshot_lake_shoreline_rank" in z.files
        assert "field__selected_snapshot_delta_fan_rank" in z.files
        assert "field__selected_snapshot_drainage_basin_id" in z.files
        assert "field__selected_snapshot_river_rank" in z.files
        assert "field__selected_snapshot_river_path_rank" in z.files
        assert "mask__selected_snapshot_lakes" in z.files
        assert "mask__selected_snapshot_deltas" in z.files
        assert "mask__selected_snapshot_delta_plain" in z.files
        assert "mask__selected_snapshot_reef_atoll" in z.files
        assert "mask__selected_snapshot_marine_shoal" in z.files
        assert "mask__selected_snapshot_seamount_shoal" in z.files
        assert "mask__selected_snapshot_oceanic_plateau_shoal" in z.files
        assert "mask__selected_snapshot_microcontinent_shoal" in z.files
        assert "mask__selected_snapshot_island_arc_shoal" in z.files
        assert "mask__selected_snapshot_inland_seaway_tuned" in z.files
        assert "mask__selected_snapshot_process_island_candidate" in z.files
        assert "mask__selected_snapshot_atoll_candidate" in z.files
        assert "mask__selected_snapshot_process_island_promoted" in z.files
        assert "mask__selected_snapshot_atoll_islet_promoted" in z.files
        assert z["field__selected_snapshot_marine_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_shelf_slope_microrelief_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_deep_ocean_fabric_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_fracture_zone_rank"].shape == (540,)
        assert z["field__selected_snapshot_abyssal_plain_fabric_rank"].shape == (540,)
        assert z["field__selected_snapshot_island_candidate_rank"].shape == (540,)
        assert z["field__selected_snapshot_atoll_candidate_rank"].shape == (540,)
        assert z["field__selected_snapshot_reef_rim_rank"].shape == (540,)
        assert z["field__selected_snapshot_atoll_lagoon_rank"].shape == (540,)
        assert z["field__selected_snapshot_fringing_reef_rank"].shape == (540,)
        assert z["field__selected_snapshot_process_island_promotion_rank"].shape == (540,)
        assert z["field__selected_snapshot_process_island_promotion_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_islet_microshape_rank"].shape == (540,)
        assert z["field__selected_snapshot_atoll_microshape_rank"].shape == (540,)
        assert z["field__selected_snapshot_island_atoll_microrelief_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_submarine_highland_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_inland_seaway_tuning_delta_m"].shape == (540,)
        assert z["field__selected_snapshot_inland_seaway_tuning_rank"].shape == (540,)
        assert z["mask__selected_snapshot_inland_seaway_landback"].shape == (540,)
        assert z["field__selected_snapshot_seamount_peak_rank"].shape == (540,)
        assert z["field__selected_snapshot_seamount_apron_rank"].shape == (540,)
        assert z["field__selected_snapshot_oceanic_plateau_edge_rank"].shape == (540,)
        assert z["field__selected_snapshot_abyssal_hill_field_rank"].shape == (540,)
        assert z["field__selected_snapshot_fluvial_microrelief_delta_m"].shape == (540,)
        assert z[
            "field__selected_snapshot_lowland_alluvial_microrelief_delta_m"
        ].shape == (540,)
        assert z["field__selected_snapshot_alluvial_fan_rank"].shape == (540,)
        assert z["field__selected_snapshot_lowland_plain_rank"].shape == (540,)
        assert z["field__selected_snapshot_piedmont_apron_rank"].shape == (540,)
        assert z["field__selected_snapshot_coastal_delta_m"].shape == (540,)
        assert z[
            "field__selected_snapshot_coastal_depositional_microrelief_delta_m"
        ].shape == (540,)
        assert z["field__selected_snapshot_coastal_plain_rank"].shape == (540,)
        assert z["field__selected_snapshot_coastal_cliff_rank"].shape == (540,)
        assert z["field__selected_snapshot_shoreface_rank"].shape == (540,)
        assert z["field__selected_snapshot_barrier_lagoon_rank"].shape == (540,)
        assert z["field__selected_snapshot_estuary_rank"].shape == (540,)
        assert z["field__selected_snapshot_delta_distributary_rank"].shape == (540,)
        assert z["field__selected_snapshot_estuary_funnel_rank"].shape == (540,)
        assert z["field__selected_snapshot_barrier_spit_rank"].shape == (540,)
        assert z["field__selected_snapshot_delta_mouth_bar_rank"].shape == (540,)
        assert z["field__selected_snapshot_estuary_tidal_channel_rank"].shape == (540,)
        assert z["field__selected_snapshot_coastal_depositional_plain_rank"].shape == (540,)
        assert z["field__selected_snapshot_strandplain_rank"].shape == (540,)
        assert z["field__selected_snapshot_tidal_flat_rank"].shape == (540,)
        assert z["field__selected_snapshot_basin_trunk_rank"].shape == (540,)
        assert z["field__selected_snapshot_floodplain_rank"].shape == (540,)
        assert z["field__selected_snapshot_meander_belt_rank"].shape == (540,)
        assert z["field__selected_snapshot_meander_scroll_rank"].shape == (540,)
        assert z["field__selected_snapshot_floodplain_swale_rank"].shape == (540,)
        assert z["field__selected_snapshot_lake_basin_rank"].shape == (540,)
        assert z["field__selected_snapshot_lake_shoreline_rank"].shape == (540,)
        assert z["field__selected_snapshot_delta_fan_rank"].shape == (540,)
        assert z["field__selected_snapshot_drainage_basin_id"].shape == (540,)
        assert z["field__selected_snapshot_river_rank"].shape == (540,)
        assert z["field__selected_snapshot_river_path_rank"].shape == (540,)
        assert z["field__selected_snapshot_fluvial_microrelief_delta_m"].shape == (540,)
        lake_core = z["mask__selected_snapshot_lakes"].astype(bool)
        lake_shore = z["field__selected_snapshot_lake_shoreline_rank"]
        if np.any(lake_core):
            assert np.max(lake_shore[lake_core]) == 0.0
        assert z[
            "field__selected_snapshot_coastal_process_microrelief_delta_m"
        ].shape == (540,)
        delta_distributary = z["field__selected_snapshot_delta_distributary_rank"]
        coastal_process = z[
            "field__selected_snapshot_coastal_process_microrelief_delta_m"
        ]
        delta_cells = delta_distributary > 0.20
        if np.any(delta_cells):
            assert float(np.nanmean(coastal_process[delta_cells])) >= 0.0

    saved = json.loads(
        (tmp_path / "refined" / "selected_snapshot_refinement_metrics.json").read_text()
    )
    assert saved["assets"]["refinement"]["refinement_delta_m.png"].endswith(
        "refinement_delta_m.png"
    )
    assert saved["assets"]["refinement"]["selected_snapshot_hydrology.png"].endswith(
        "selected_snapshot_hydrology.png"
    )
    assert saved["assets"]["refinement"]["selected_snapshot_drainage_basins.png"].endswith(
        "selected_snapshot_drainage_basins.png"
    )
    assert saved["assets"]["refinement"][
        "selected_snapshot_fluvial_lacustrine_microshapes.png"
    ].endswith("selected_snapshot_fluvial_lacustrine_microshapes.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png"
    ].endswith("selected_snapshot_fluvial_lacustrine_microshapes_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_fluvial_microrelief_delta.png"
    ].endswith("selected_snapshot_fluvial_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_fluvial_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_fluvial_microrelief_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_lowland_alluvial_microrelief_delta.png"
    ].endswith("selected_snapshot_lowland_alluvial_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_lowland_alluvial_microrelief_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_marine_microgeomorphology.png"
    ].endswith("selected_snapshot_marine_microgeomorphology.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_marine_object_classes.png"
    ].endswith("selected_snapshot_marine_object_classes.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_shelf_slope_microrelief_delta.png"
    ].endswith("selected_snapshot_shelf_slope_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_shelf_slope_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_shelf_slope_microrelief_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_deep_ocean_fabric_delta.png"
    ].endswith("selected_snapshot_deep_ocean_fabric_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_deep_ocean_fabric_zoom_sheet.png"
    ].endswith("selected_snapshot_deep_ocean_fabric_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_submarine_highlands.png"
    ].endswith("selected_snapshot_submarine_highlands.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_submarine_highlands_zoom_sheet.png"
    ].endswith("selected_snapshot_submarine_highlands_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_island_atoll_candidates.png"
    ].endswith("selected_snapshot_island_atoll_candidates.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_process_island_promotion.png"
    ].endswith("selected_snapshot_process_island_promotion.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_process_island_promotion_zoom_sheet.png"
    ].endswith("selected_snapshot_process_island_promotion_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_island_atoll_microshapes.png"
    ].endswith("selected_snapshot_island_atoll_microshapes.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_island_atoll_microshapes_zoom_sheet.png"
    ].endswith("selected_snapshot_island_atoll_microshapes_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_island_atoll_microrelief_delta.png"
    ].endswith("selected_snapshot_island_atoll_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_island_atoll_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_island_atoll_microrelief_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_reef_atoll_morphology.png"
    ].endswith("selected_snapshot_reef_atoll_morphology.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_reef_atoll_zoom_sheet.png"
    ].endswith("selected_snapshot_reef_atoll_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_morphology.png"
    ].endswith("selected_snapshot_coastal_morphology.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_zoom_sheet.png"
    ].endswith("selected_snapshot_coastal_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_process_linework.png"
    ].endswith("selected_snapshot_coastal_process_linework.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_process_linework_zoom_sheet.png"
    ].endswith("selected_snapshot_coastal_process_linework_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_process_microrelief_delta.png"
    ].endswith("selected_snapshot_coastal_process_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_process_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_coastal_process_microrelief_zoom_sheet.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_depositional_microrelief_delta.png"
    ].endswith("selected_snapshot_coastal_depositional_microrelief_delta.png")
    assert saved["assets"]["refinement"][
        "selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png"
    ].endswith("selected_snapshot_coastal_depositional_microrelief_zoom_sheet.png")


def test_selected_snapshot_process_island_promotion_is_explicit(tmp_path):
    grid = SphereGrid.fibonacci(240, CONSTANTS.EARTH_RADIUS)
    land = (grid.lat > 34.0) | ((grid.lat < -30.0) & (grid.lon < -110.0))
    rel = np.where(land, 220.0 + 3.0 * np.abs(grid.lat), -3200.0)
    process_high = (
        ~land
        & (np.abs(grid.lat) < 18.0)
        & (grid.lon > -70.0)
        & (grid.lon < 40.0)
    )
    rel[process_high] = -90.0
    shallow_rim = (
        ~land
        & (np.abs(grid.lat - 8.0) < 8.0)
        & (grid.lon > 65.0)
        & (grid.lon < 110.0)
    )
    rel[shallow_rim] = -70.0
    detail = np.where(land, CONT_DETAIL_PLATFORM, 0).astype(np.int32)
    depth_province = np.where(land, 0, OCEAN_DEPTH_ABYSS).astype(np.int32)
    arrays_path = tmp_path / "p107_terminal_arrays.npz"
    np.savez_compressed(
        arrays_path,
        grid_lat=grid.lat.astype(np.float32),
        grid_lon=grid.lon.astype(np.float32),
        grid_cell_area_m2=grid.cell_area.astype(np.float64),
        field__terrain_elevation_m=rel.astype(np.float32),
        field__tectonics_plate_id=(np.arange(grid.n) % 6).astype(np.int32),
        field__crust_type=land.astype(np.int32),
        field__crust_age_myr=np.where(land, 1500.0, 45.0).astype(np.float32),
        field__ocean_depth_province=depth_province,
        field__terrain_continental_detail=detail,
        field__terrain_continental_detail_region_code=detail,
        field__terrain_orogenic_parent_hierarchy=np.zeros(grid.n, dtype=np.int32),
        field__terrain_orogenic_hierarchy_spine=np.zeros(grid.n, dtype=np.int32),
        field__terrain_orogenic_shoulder_halo=np.zeros(grid.n, dtype=np.int32),
        field__terrain_orogenic_highland_apron=np.zeros(grid.n, dtype=np.int32),
        boundary__ridge=np.zeros(grid.n, dtype=np.uint8),
        boundary__trench=np.zeros(grid.n, dtype=np.uint8),
        object__province_mid_ocean_ridge=np.zeros(grid.n, dtype=np.uint8),
        object__margin_trench=np.zeros(grid.n, dtype=np.uint8),
        object__seamount_chain=process_high.astype(np.uint8),
        object__oceanic_plateau=shallow_rim.astype(np.uint8),
        object__microcontinent=(process_high | shallow_rim).astype(np.uint8),
        object__island_arc=shallow_rim.astype(np.uint8),
        object__mountain_orogen=np.zeros(grid.n, dtype=np.uint8),
        object__parent_orogen_crest=np.zeros(grid.n, dtype=np.uint8),
    )
    metrics = {
        "schema": "aevum.p107_terminal_world_audit.v1",
        "cells": grid.n,
        "sea_level_m": 0.0,
        "array_archive": {
            "path": str(arrays_path),
            "manifest": {
                "fields": {
                    "terrain.elevation_m": "field__terrain_elevation_m",
                    "tectonics.plate_id": "field__tectonics_plate_id",
                    "crust.type": "field__crust_type",
                    "crust.age_myr": "field__crust_age_myr",
                    "ocean.depth_province": "field__ocean_depth_province",
                    "terrain.continental_detail": "field__terrain_continental_detail",
                    "terrain.continental_detail_region_code": (
                        "field__terrain_continental_detail_region_code"
                    ),
                    "terrain.orogenic_parent_hierarchy": (
                        "field__terrain_orogenic_parent_hierarchy"
                    ),
                    "terrain.orogenic_hierarchy_spine": (
                        "field__terrain_orogenic_hierarchy_spine"
                    ),
                    "terrain.orogenic_shoulder_halo": (
                        "field__terrain_orogenic_shoulder_halo"
                    ),
                    "terrain.orogenic_highland_apron": (
                        "field__terrain_orogenic_highland_apron"
                    ),
                }
            },
        },
        "p109_hypsometry_comparison": {"within_p109_envelope": True},
        "p110a_modern_planform": {
            "summary": {"land_fraction": float(np.mean(land))},
            "warning_flags": [],
        },
    }
    metrics_path = tmp_path / "p107_terminal_metrics.json"
    metrics_path.write_text(json.dumps(metrics) + "\n")

    default_summary = refine_selected_snapshot(
        SelectedSnapshotRefinementConfig(
            source=metrics_path,
            outdir=tmp_path / "default_refined",
            target_cells=720,
            width=72,
            height=36,
            detail_seed=2026,
        )
    )
    promoted_summary = refine_selected_snapshot(
        SelectedSnapshotRefinementConfig(
            source=metrics_path,
            outdir=tmp_path / "promoted_refined",
            target_cells=720,
            width=72,
            height=36,
            detail_seed=2026,
            allow_process_islands=True,
            render_groups=("island-atoll",),
        )
    )

    assert default_summary["land_ocean_sign_flip_fraction"] == 0.0
    assert default_summary["process_island_promoted_cell_fraction_parent_ocean"] == 0.0
    assert promoted_summary["allow_process_islands"] is True
    assert promoted_summary["render_groups"] == ["island-atoll"]
    assert promoted_summary["land_ocean_sign_flip_fraction"] > 0.0
    assert promoted_summary["process_island_promoted_cell_fraction_parent_ocean"] > 0.0
    assert promoted_summary["process_island_promotion_delta_nonzero_abs_p95_m"] > 0.0
    refinement_assets = promoted_summary["assets"]["refinement"]
    assert "selected_snapshot_island_atoll_microrelief_zoom_sheet.png" in refinement_assets
    assert "selected_snapshot_reef_atoll_zoom_sheet.png" in refinement_assets
    assert "selected_snapshot_hydrology.png" not in refinement_assets
    assert promoted_summary["assets"]["p107"] == {}
    render_summary = render_selected_snapshot_refinement_assets(
        tmp_path / "promoted_refined",
        render_groups=("submarine",),
        width=72,
        height=36,
        outdir=tmp_path / "promoted_submarine_rendered",
    )
    assert render_summary["render_groups"] == ["submarine"]
    assert (
        "selected_snapshot_submarine_highlands_zoom_sheet.png"
        in render_summary["assets"]["refinement"]
    )
    assert "selected_snapshot_hydrology.png" not in render_summary["assets"]["refinement"]
    with np.load(
        tmp_path / "promoted_refined" / "selected_snapshot_refined_arrays.npz"
    ) as z:
        promoted = z["mask__selected_snapshot_process_island_promoted"].astype(bool)
        refined = z["field__terrain_elevation_m"]
        parent = z["field__selected_snapshot_parent_elevation_rel_m"]
        islet_microshape = z["field__selected_snapshot_islet_microshape_rank"]
        atoll_microshape = z["field__selected_snapshot_atoll_microshape_rank"]
        assert np.count_nonzero(promoted) > 0
        assert np.all(parent[promoted] < 0.0)
        assert np.all(refined[promoted] > 0.0)
        assert np.any(islet_microshape[promoted] > 0.0)
        assert np.count_nonzero(atoll_microshape > 0.0) >= 0
