"""Visualisation helpers (matplotlib, Agg backend).

Renders truth-layer fields, the compiled hex map, the deep-time event timeline,
and an elevation history strip.  Rendering is deliberately the LAST thing built;
it never feeds back into the simulation.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import (
    BoundaryNorm,
    FuncNorm,
    LinearSegmentedColormap,
    ListedColormap,
    TwoSlopeNorm,
)

from aevum.compiler.map_compiler import CompiledMap, TERRAIN_NAMES
from aevum.modules.tectonics import INTERNAL_BLOCK_NAMES

BIOME_COLORS = [
    "#173b57",  # ocean
    "#eef5f8",  # ice
    "#d9c98f",  # desert
    "#a8b96f",  # grassland
    "#3f7f4a",  # forest
    "#b9c6b8",  # tundra
    "#176c43",  # tropical
]
TERRAIN_COLORS = [
    "#173b57",  # ocean
    "#5f9fb1",  # coast / shallow sea
    "#3f8a45",  # plains
    "#d88a32",  # hills
    "#6f3f8f",  # mountains
    "#eef5f8",  # ice
]
ELEVATION_DISPLAY_MIN_M = -6000.0
ELEVATION_DISPLAY_MAX_M = 6000.0
ELEVATION_COLOR_TICKS = [
    -6000, -4500, -3000, -1500, 0, 200, 500, 1000, 2000, 3000, 4500, 6000
]
ELEVATION_COLOR_VALUES = np.array(
    [-6000.0, -4500.0, -3000.0, -1500.0, -1.0, 0.0,
     200.0, 500.0, 1000.0, 2000.0, 3000.0, 4500.0, 5200.0, 6000.0],
    dtype=np.float64,
)
ELEVATION_COLOR_POSITIONS = np.array(
    [0.00, 0.10, 0.21, 0.33, 0.44, 0.45,
     0.53, 0.61, 0.69, 0.77, 0.85, 0.93, 0.965, 1.00],
    dtype=np.float64,
)


def _elevation_forward(values):
    return np.interp(
        values,
        ELEVATION_COLOR_VALUES,
        ELEVATION_COLOR_POSITIONS,
        left=0.0,
        right=1.0,
    )


def _elevation_inverse(values):
    return np.interp(
        values,
        ELEVATION_COLOR_POSITIONS,
        ELEVATION_COLOR_VALUES,
        left=ELEVATION_DISPLAY_MIN_M,
        right=ELEVATION_DISPLAY_MAX_M,
    )


def _elevation_color_stop(elevation_m: float) -> float:
    return float(_elevation_forward(float(elevation_m)))


ELEVATION_CMAP = LinearSegmentedColormap.from_list(
    "aevum_elevation_absolute",
    [
        (_elevation_color_stop(-6000.0), "#06162e"),
        (_elevation_color_stop(-4500.0), "#0b355e"),
        (_elevation_color_stop(-3000.0), "#165f7b"),
        (_elevation_color_stop(-1500.0), "#4b9aa7"),
        (_elevation_color_stop(-1.0), "#b7e3dc"),
        (_elevation_color_stop(0.0), "#064b2c"),
        (_elevation_color_stop(200.0), "#1f6e3c"),
        (_elevation_color_stop(500.0), "#78ad55"),
        (_elevation_color_stop(1000.0), "#c4d85c"),
        (_elevation_color_stop(2000.0), "#f0dd63"),
        (_elevation_color_stop(3000.0), "#df8c34"),
        (_elevation_color_stop(4500.0), "#c43b32"),
        (_elevation_color_stop(5200.0), "#b58cc8"),
        (_elevation_color_stop(6000.0), "#fbfbf4"),
    ],
)
ELEVATION_NORM = FuncNorm(
    (_elevation_forward, _elevation_inverse),
    vmin=ELEVATION_DISPLAY_MIN_M,
    vmax=ELEVATION_DISPLAY_MAX_M,
    clip=True,
)
HEX_LAND_ELEVATION_COLORS = [
    ("lowland", "#1f5f38"),
    ("plains", "#5f9f4a"),
    ("upland", "#a9c94f"),
    ("high plain", "#efe89a"),
    ("hills", "#d98a34"),
    ("highlands", "#bd3f32"),
    ("mountains", "#6f3f8f"),
    ("alpine", "#f5f5ef"),
]
TEMPERATURE_CMAP = LinearSegmentedColormap.from_list(
    "aevum_temperature",
    ["#24476b", "#6e9fbe", "#e6ece8", "#d89b73", "#8f2f2f"],
)
PRECIP_CMAP = LinearSegmentedColormap.from_list(
    "aevum_precip",
    ["#efe9c7", "#b7ce8a", "#78b6a6", "#3f8fb7", "#1c4f86"],
)
OCEAN_CRUST_AGE_CMAP = LinearSegmentedColormap.from_list(
    "aevum_ocean_crust_age",
    ["#12243f", "#1f5d7a", "#4fa5a0", "#bfd07b", "#d79a5b"],
)
CONTINENT_CRUST_AGE_CMAP = LinearSegmentedColormap.from_list(
    "aevum_continent_crust_age",
    ["#ded8c8", "#c7ab82", "#997160", "#6f5149"],
)


def to_raster(grid, values, width=360, height=180):
    lon = np.linspace(-180, 180, width)
    lat = np.linspace(90, -90, height)
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    idx = grid.nearest_latlon(LAT.ravel(), LON.ravel()).reshape(height, width)
    return values[idx]


def to_raster_continuous(
    grid,
    values,
    width=360,
    height=180,
    *,
    k: int = 6,
    preserve_sign: bool = False,
):
    """Rasterize a continuous spherical field without exposing cell rows.

    Discrete layers such as plates, biomes, and object masks must keep using
    :func:`to_raster`.  Elevation and bathymetry are continuous fields, though;
    nearest-cell projection at display resolutions much higher than the
    simulation grid makes Fibonacci cell rows look like straight geomorphic
    bands.  A small inverse-distance stencil removes that display artifact while
    optional sign preservation keeps coastlines from being blended across sea
    level.
    """
    values = np.asarray(values)
    if int(k) <= 1 or grid.n <= 1:
        return to_raster(grid, values, width=width, height=height)
    lon = np.linspace(-180.0, 180.0, int(width))
    lat = np.linspace(90.0, -90.0, int(height))
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")
    lat_rad = np.radians(lat_grid.ravel())
    lon_rad = np.radians(lon_grid.ravel())
    query = np.column_stack([
        np.cos(lat_rad) * np.cos(lon_rad),
        np.cos(lat_rad) * np.sin(lon_rad),
        np.sin(lat_rad),
    ])
    kk = max(1, min(int(k), int(grid.n)))
    dist, idx = grid._kdtree.query(query, k=kk)
    if kk == 1:
        idx = idx[:, None]
        dist = dist[:, None]
    vals = values[idx]
    nearest = vals[:, 0]
    finite = np.isfinite(vals)
    weights = 1.0 / np.maximum(dist, 1.0e-9) ** 2
    weights = np.where(finite, weights, 0.0)
    if preserve_sign:
        same_sign = np.signbit(vals) == np.signbit(nearest)[:, None]
        weights = np.where(same_sign, weights, 0.0)
    weight_sum = weights.sum(axis=1)
    out = (vals * weights).sum(axis=1) / np.maximum(weight_sum, 1.0e-12)
    bad = weight_sum <= 0.0
    out[bad] = nearest[bad]
    if preserve_sign:
        land_flip = (nearest >= 0.0) & (out < 0.0)
        ocean_flip = (nearest < 0.0) & (out >= 0.0)
        out[land_flip | ocean_flip] = nearest[land_flip | ocean_flip]
    return out.reshape(int(height), int(width))


def render_elevation_raster(ax, raster: np.ndarray, *, title: str | None = None):
    im = ax.imshow(
        raster,
        cmap=ELEVATION_CMAP,
        norm=ELEVATION_NORM,
        extent=[-180, 180, -90, 90],
    )
    if title:
        ax.set_title(title)
    return im


def add_elevation_colorbar(fig, ax, im):
    cb = fig.colorbar(im, ax=ax, shrink=0.7, ticks=ELEVATION_COLOR_TICKS)
    cb.set_label("m rel. sea level")
    return cb


def render_world(world, outdir: Path) -> list[Path]:
    grid = world.grid
    out = []
    sea = world.sea_level

    elev = world.get_field("terrain.elevation_m")
    er = to_raster_continuous(grid, elev - sea, preserve_sign=True)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = render_elevation_raster(
        ax, er, title=f"Elevation (m rel. sea level)  t={world.time_myr:.0f} Myr"
    )
    add_elevation_colorbar(fig, ax, im)
    p = outdir / "elevation.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig); out.append(p)

    for name, field, cmap, title in [
        ("temperature", world.get_field("climate.surface_temperature", 288.0) - 273.15,
         TEMPERATURE_CMAP, "Surface temperature (C)"),
        ("precip", world.get_field("climate.precipitation", 0.0), PRECIP_CMAP,
         "Precipitation (mm/yr)"),
    ]:
        r = to_raster(grid, field)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        im = ax.imshow(r, cmap=cmap, extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.7)
        p = outdir / f"{name}.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig); out.append(p)

    seasonal_temp = world.fields.get("climate.seasonal_temperature")
    if seasonal_temp is not None:
        seasonal_temp = np.asarray(seasonal_temp)
        if seasonal_temp.shape == (4, grid.n):
            labels = ("DJF", "MAM", "JJA", "SON")
            rasters = [to_raster(grid, seasonal_temp[i] - 273.15) for i in range(4)]
            combined = np.concatenate([r.ravel() for r in rasters])
            vmin = float(np.percentile(combined, 2))
            vmax = float(np.percentile(combined, 98))
            fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
            im = None
            for ax, label, raster in zip(axes.ravel(), labels, rasters):
                im = ax.imshow(raster, cmap=TEMPERATURE_CMAP, vmin=vmin, vmax=vmax,
                               extent=[-180, 180, -90, 90])
                ax.set_title(f"{label} temperature (C)")
            if im is not None:
                fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
            p = outdir / "temperature_seasons.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    for field_name, filename, title, cmap, centered, ratio_centered, kelvin_to_celsius in [
        ("climate.seasonal_sst", "seasonal_sst.png",
         "sea-surface temperature (C)", TEMPERATURE_CMAP, False, False, True),
        ("climate.seasonal_precipitation", "precip_seasons.png",
         "precipitation (mm/yr equivalent)", PRECIP_CMAP, False, False, False),
        ("cryosphere.seasonal_sea_ice", "sea_ice_seasons.png",
         "sea-ice fraction", "Blues", False, False, False),
        ("cryosphere.seasonal_snow", "snow_seasons.png",
         "snow persistence", "Blues", False, False, False),
        ("climate.seasonal_cloud_albedo_proxy", "cloud_albedo_seasons.png",
         "cloud albedo proxy", "Greys", False, False, False),
        ("atmosphere.seasonal_pressure_proxy", "seasonal_pressure.png",
         "seasonal pressure proxy", "coolwarm", True, False, False),
        ("atmosphere.moisture_access", "moisture_access.png",
         "moisture access", PRECIP_CMAP, False, False, False),
        ("atmosphere.moisture_source_basin_id", "moisture_source_basin_id.png",
         "moisture source basin id", "tab20", False, False, False),
        ("atmosphere.monsoon_potential", "monsoon_potential.png",
         "monsoon potential", "coolwarm", True, False, False),
        ("atmosphere.source_ocean_warmth", "source_ocean_warmth.png",
         "source ocean warmth", TEMPERATURE_CMAP, False, False, False),
        ("climate.monsoon_rainfall_corridor", "monsoon_rainfall_corridor.png",
         "monsoon rainfall corridor", PRECIP_CMAP, False, False, False),
        ("climate.storm_track_rainfall_corridor",
         "storm_track_rainfall_corridor.png",
         "storm-track rainfall corridor", PRECIP_CMAP, False, False, False),
        ("climate.rain_shadow_index", "rain_shadow_index.png",
         "rain-shadow index", PRECIP_CMAP, False, False, False),
        ("climate.regional_precipitation_response",
         "regional_precipitation_response.png",
         "regional precipitation response", "coolwarm", False, True, False),
        ("climate.moisture_flow_precipitation_response",
         "moisture_flow_precipitation_response.png",
         "moisture-flow precipitation response", "coolwarm", False, True, False),
        ("climate.moisture_budget_region_id",
         "moisture_budget_region_id.png",
         "moisture budget region id", "tab20", False, False, False),
        ("climate.precipitation_response_region_id",
         "precipitation_response_region_id.png",
         "precipitation response region id", "tab20", False, False, False),
        ("climate.receiver_catchment_id",
         "receiver_catchment_id.png",
         "receiver catchment id", "tab20", False, False, False),
        ("climate.source_basin_supply_index",
         "source_basin_supply_index.png",
         "source-basin supply index", PRECIP_CMAP, False, False, False),
        ("climate.receiver_catchment_supply_balance",
         "receiver_catchment_supply_balance.png",
         "receiver catchment supply balance", PRECIP_CMAP, False, False, False),
        ("climate.receiver_supply_precipitation_feedback",
         "receiver_supply_precipitation_feedback.png",
         "receiver-supply precipitation feedback", "coolwarm", False, True, False),
    ]:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field)
        if field.shape != (4, grid.n):
            continue
        if kelvin_to_celsius:
            field = field - 273.15
        labels = ("DJF", "MAM", "JJA", "SON")
        rasters = [to_raster(grid, field[i]) for i in range(4)]
        combined = np.concatenate([r.ravel() for r in rasters])
        fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
        im = None
        norm = None
        if ratio_centered:
            span = max(float(np.nanpercentile(np.abs(combined - 1.0), 98)), 0.05)
            vmin = max(0.0, 1.0 - span)
            vmax = 1.0 + span
            norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)
        elif centered:
            vmax = max(float(np.percentile(np.abs(combined), 98)), 0.1)
            vmin = -vmax
        else:
            vmin = min(0.0, float(np.percentile(combined, 2)))
            vmax = max(float(np.percentile(combined, 98)), 0.1)
        for ax, label, raster in zip(axes.ravel(), labels, rasters):
            im = ax.imshow(
                raster,
                cmap=cmap,
                norm=norm,
                vmin=None if norm is not None else vmin,
                vmax=None if norm is not None else vmax,
                extent=[-180, 180, -90, 90],
            )
            ax.set_title(f"{label} {title}")
        if im is not None:
            fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
        p = outdir / filename
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        out.append(p)

    seasonal_wind = world.fields.get("atmosphere.seasonal_wind")
    if seasonal_wind is not None:
        seasonal_wind = np.asarray(seasonal_wind)
        if seasonal_wind.shape == (4, grid.n, 3):
            labels = ("DJF", "MAM", "JJA", "SON")
            speed = np.linalg.norm(seasonal_wind, axis=2)
            rasters = [to_raster(grid, speed[i]) for i in range(4)]
            combined = np.concatenate([r.ravel() for r in rasters])
            vmax = float(np.percentile(combined, 98)) if combined.size else 1.0
            fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
            im = None
            for ax, label, raster in zip(axes.ravel(), labels, rasters):
                im = ax.imshow(raster, cmap="viridis", vmin=0.0, vmax=max(vmax, 1.0),
                               extent=[-180, 180, -90, 90])
                ax.set_title(f"{label} wind speed (m/s)")
            if im is not None:
                fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
            p = outdir / "wind_seasons.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    currents = world.fields.get("ocean.currents")
    if currents is not None:
        currents = np.asarray(currents)
        if currents.shape == (grid.n, 3):
            speed = np.linalg.norm(currents, axis=1)
            r = to_raster(grid, speed)
            vmax = float(np.percentile(speed[speed > 0.0], 98)) if (speed > 0.0).any() else 1.0
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(r, cmap="cividis", vmin=0.0, vmax=max(vmax, 0.1),
                           extent=[-180, 180, -90, 90])
            ax.set_title("Ocean current speed (m/s)")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "currents.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    pressure_proxy = world.fields.get("atmosphere.land_sea_pressure_proxy")
    if pressure_proxy is not None:
        pressure_proxy = np.asarray(pressure_proxy)
        if pressure_proxy.shape == (4, grid.n):
            labels = ("DJF", "MAM", "JJA", "SON")
            rasters = [to_raster(grid, pressure_proxy[i]) for i in range(4)]
            combined = np.concatenate([r.ravel() for r in rasters])
            vmax = max(float(np.percentile(np.abs(combined), 98)), 0.1)
            fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
            im = None
            for ax, label, raster in zip(axes.ravel(), labels, rasters):
                im = ax.imshow(raster, cmap="coolwarm", vmin=-vmax, vmax=vmax,
                               extent=[-180, 180, -90, 90])
                ax.set_title(f"{label} land-sea pressure proxy")
            if im is not None:
                fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
            p = outdir / "land_sea_pressure.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    for field_name, filename, title in [
        ("atmosphere.thermal_wind_anomaly", "thermal_wind_anomaly.png",
         "thermal wind anomaly"),
        ("atmosphere.orographic_wind_anomaly", "orographic_wind_anomaly.png",
         "orographic wind anomaly"),
    ]:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field)
        if field.shape != (4, grid.n, 3):
            continue
        labels = ("DJF", "MAM", "JJA", "SON")
        speed = np.linalg.norm(field, axis=2)
        rasters = [to_raster(grid, speed[i]) for i in range(4)]
        combined = np.concatenate([r.ravel() for r in rasters])
        vmax = float(np.percentile(combined, 98)) if combined.size else 1.0
        fig, axes = plt.subplots(2, 2, figsize=(10, 5.6), constrained_layout=True)
        im = None
        for ax, label, raster in zip(axes.ravel(), labels, rasters):
            im = ax.imshow(raster, cmap="magma", vmin=0.0, vmax=max(vmax, 0.1),
                           extent=[-180, 180, -90, 90])
            ax.set_title(f"{label} {title} (m/s)")
        if im is not None:
            fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
        p = outdir / filename
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        out.append(p)

    itcz_lat = world.fields.get("atmosphere.itcz_latitude")
    if itcz_lat is not None:
        itcz_lat = np.asarray(itcz_lat)
        if itcz_lat.shape == (4,):
            labels = ("DJF", "MAM", "JJA", "SON")
            colors = ("#2c7fb8", "#41ab5d", "#d95f0e", "#756bb1")
            fig, ax = plt.subplots(figsize=(9, 4.5))
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            ax.axhline(0.0, color="#666666", linewidth=0.8, alpha=0.6)
            for lat, label, color in zip(itcz_lat, labels, colors):
                ax.axhline(float(lat), color=color, linewidth=2.0, label=label)
            ax.set_title("Seasonal ITCZ latitude")
            ax.set_xlabel("longitude")
            ax.set_ylabel("latitude")
            ax.legend(loc="upper right", frameon=False)
            p = outdir / "itcz_track.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    optional_climate_maps = [
        ("temperature_seasonality", "climate.temperature_seasonality",
         TEMPERATURE_CMAP, "Temperature seasonality (C)"),
        ("precip_seasonality", "climate.precipitation_seasonality",
         PRECIP_CMAP, "Precipitation seasonality"),
        ("sea_ice", "cryosphere.sea_ice", "Blues", "Sea-ice fraction"),
        ("snow_ice", "cryosphere.snow_persistence", "Blues",
         "Snow persistence"),
        ("cloud_albedo", "climate.cloud_albedo_proxy", "Greys",
         "Cloud albedo proxy"),
        ("vegetation_climate_feedback", "biosphere.vegetation_climate_feedback",
         "YlGn", "Vegetation climate feedback proxy"),
        ("monsoon_index", "climate.monsoon_index", PRECIP_CMAP, "Monsoon index"),
        ("dry_season_length", "climate.dry_season_length", PRECIP_CMAP,
         "Dry season length"),
        ("wet_season_peak", "climate.wet_season_peak", "tab10", "Wet season peak"),
        ("moisture_convergence", "climate.moisture_convergence", PRECIP_CMAP,
         "Moisture convergence"),
        ("evaporation", "climate.evaporation", PRECIP_CMAP, "Evaporation"),
        ("runoff", "climate.runoff", PRECIP_CMAP, "Runoff"),
        ("orographic_precipitation", "climate.orographic_precipitation",
         PRECIP_CMAP, "Orographic precipitation"),
        ("continent_id", "climate.continent_id", "tab20", "Continent id"),
        ("continent_interiority", "climate.continent_interiority",
         PRECIP_CMAP, "Continent interiority"),
        ("coast_distance", "climate.coast_distance", PRECIP_CMAP,
         "Distance from coast"),
        ("coast_strength", "climate.coast_strength", PRECIP_CMAP,
         "Coast strength"),
        ("coast_orientation", "climate.coast_facing_east", "coolwarm",
         "Coast orientation east-west"),
        ("ocean_heat_flux", "climate.ocean_heat_flux",
         TEMPERATURE_CMAP, "Ocean heat flux"),
        ("coupling_residual", "climate.coupling_residual",
         PRECIP_CMAP, "SST/wind coupling residual"),
        ("ocean_heat_transport", "ocean.current_heat_transport",
         TEMPERATURE_CMAP, "Ocean heat transport proxy"),
        ("sst_anomaly", "ocean.sst_anomaly", TEMPERATURE_CMAP,
         "Ocean SST anomaly"),
        ("current_streamfunction", "ocean.current_streamfunction", "coolwarm",
         "Ocean current streamfunction"),
        ("gyres", "ocean.gyre_id", "tab20", "Ocean gyre id"),
        ("boundary_current_type", "ocean.boundary_current_type", "coolwarm",
         "Boundary current type"),
        ("strait_exchange", "ocean.strait_exchange", PRECIP_CMAP,
         "Strait exchange"),
        ("upwelling", "ocean.upwelling", PRECIP_CMAP, "Coastal upwelling proxy"),
        ("ocean_basin_id", "ocean.basin_id", "tab20", "Ocean basin id"),
        ("ocean_depth_province", "ocean.depth_province", "tab20",
         "Ocean depth province"),
        ("ocean_margin_type", "ocean.margin_type", "tab20", "Ocean margin type"),
        ("ocean_shelf_width", "ocean.shelf_width", PRECIP_CMAP, "Ocean shelf width"),
        ("ocean_gateway_id", "ocean.gateway_id", "tab20", "Ocean gateway id"),
        ("shelf_index", "ocean.shelf_index", PRECIP_CMAP, "Ocean shelf index"),
        ("strait_index", "ocean.strait_index", PRECIP_CMAP, "Ocean strait index"),
        ("terrain_barriers", "terrain.barrier_index", PRECIP_CMAP,
         "Terrain barrier index"),
        ("wind_gaps", "terrain.wind_gap_index", PRECIP_CMAP, "Wind gap index"),
        ("terrain_blocking", "atmosphere.terrain_blocking", PRECIP_CMAP,
         "Terrain blocking"),
        ("geographic_circulation_index", "atmosphere.geographic_circulation_index",
         PRECIP_CMAP, "Geographic circulation index"),
    ]
    for name, field_name, cmap, title in optional_climate_maps:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field)
        if field.shape != (grid.n,):
            continue
        r = to_raster(grid, field)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        if name in {"ocean_heat_flux", "ocean_heat_transport", "sst_anomaly"}:
            vmax = max(float(np.percentile(np.abs(field), 98)), 0.1)
            im = ax.imshow(r, cmap="coolwarm", vmin=-vmax, vmax=vmax,
                           extent=[-180, 180, -90, 90])
        elif name in {"ocean_basin_id", "continent_id"}:
            im = ax.imshow(r, cmap=cmap, vmin=-1.0, vmax=max(float(np.nanmax(field)), 1.0),
                           extent=[-180, 180, -90, 90])
        elif name == "coast_orientation":
            im = ax.imshow(r, cmap=cmap, vmin=-1.0, vmax=1.0,
                           extent=[-180, 180, -90, 90])
        elif name in {"monsoon_index"}:
            vmax = max(float(np.percentile(np.abs(field), 98)), 0.1)
            im = ax.imshow(r, cmap="coolwarm", vmin=-vmax, vmax=vmax,
                           extent=[-180, 180, -90, 90])
        elif name == "wet_season_peak":
            im = ax.imshow(r, cmap=cmap, vmin=0.0, vmax=3.0,
                           extent=[-180, 180, -90, 90])
        else:
            im = ax.imshow(r, cmap=cmap, extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.7)
        p = outdir / f"{name}.png"
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        out.append(p)

    p = _render_geography_primitives(world, outdir)
    if p is not None:
        out.append(p)

    out.extend(_render_ocean_geography(world, outdir))
    p = _render_seam_continuity(world, outdir)
    if p is not None:
        out.append(p)
    out.extend(_render_tectonic_object_layers(world, outdir))
    out.extend(_render_morphology_diagnostics(world, outdir))
    p = _render_terrain_province(world, outdir)
    if p is not None:
        out.append(p)
    p = _render_continental_detail(world, outdir)
    if p is not None:
        out.append(p)
    p = _render_inland_geomorphology_regions(world, outdir)
    if p is not None:
        out.append(p)
    out.extend(_render_internal_geographic_blocks(world, outdir))
    p = _render_wilson_cycle_phase(world, outdir)
    if p is not None:
        out.append(p)
    p = _render_ocean_gateways(world, outdir)
    if p is not None:
        out.append(p)

    p = _render_crust_age(world, outdir)
    out.append(p)

    # plates (discrete)
    plate = world.get_field("tectonics.plate_id", 0.0)
    r = to_raster(grid, plate)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.imshow(r, cmap="tab20", extent=[-180, 180, -90, 90])
    ax.set_title("Plates")
    p = outdir / "plates.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig); out.append(p)

    # biome (discrete)
    biome = world.get_field("biosphere.biome", 0.0)
    r = to_raster(grid, biome)
    cmap = ListedColormap(BIOME_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(BIOME_COLORS) + 0.5), cmap.N)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.imshow(r, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Biomes")
    p = outdir / "biomes.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig); out.append(p)
    return out


def _render_crust_age(world, outdir: Path) -> Path:
    """Render oceanic spreading age without ancient continents crushing contrast."""
    grid = world.grid
    age = world.get_field("crust.age_myr", 0.0)
    ctype = world.get_field("crust.type", 0.0)
    age_r = to_raster(grid, age)
    type_r = to_raster(grid, ctype)
    ocean = type_r < 0.5
    continent = ~ocean

    ocean_age = np.ma.masked_where(~ocean, age_r)
    cont_age = np.ma.masked_where(~continent, age_r)
    ocean_values = age[ctype < 0.5]
    cont_values = age[ctype >= 0.5]
    ocean_vmax = 300.0
    if ocean_values.size:
        ocean_vmax = float(np.clip(np.percentile(ocean_values, 98), 120.0, 500.0))
    cont_vmax = 3000.0
    if cont_values.size:
        cont_vmax = float(np.clip(np.percentile(cont_values, 98), 900.0, 3500.0))

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.set_facecolor("#efe8d7")
    ax.imshow(cont_age, cmap=CONTINENT_CRUST_AGE_CMAP, vmin=0.0, vmax=cont_vmax,
              extent=[-180, 180, -90, 90], alpha=0.58)
    im = ax.imshow(ocean_age, cmap=OCEAN_CRUST_AGE_CMAP, vmin=0.0, vmax=ocean_vmax,
                   extent=[-180, 180, -90, 90])
    ax.set_title("Oceanic crust age (Myr); continents muted")
    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label("oceanic crust age (Myr)")
    p = outdir / "crust_age.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_geography_primitives(world, outdir: Path) -> Path | None:
    grid = world.grid
    specs = [
        ("climate.continent_id", "Continents", "tab20", None),
        ("ocean.basin_id", "Ocean basins", "tab20", None),
        ("ocean.shelf_index", "Shelves", PRECIP_CMAP, (0.0, 1.0)),
        ("ocean.strait_index", "Straits", PRECIP_CMAP, (0.0, 1.0)),
        ("terrain.barrier_index", "Terrain barriers", PRECIP_CMAP, (0.0, 1.0)),
        ("terrain.wind_gap_index", "Wind gaps", PRECIP_CMAP, (0.0, 1.0)),
    ]
    if not any(field in world.fields for field, *_ in specs):
        return None

    fig, axes = plt.subplots(3, 2, figsize=(10, 8.2), constrained_layout=True)
    for ax, (field_name, title, cmap, limits) in zip(axes.ravel(), specs):
        field = world.fields.get(field_name)
        if field is None:
            ax.axis("off")
            continue
        field = np.asarray(field)
        if field.shape != (grid.n,):
            ax.axis("off")
            continue
        r = to_raster(grid, field)
        if field_name in {"climate.continent_id", "ocean.basin_id"}:
            vmax = max(float(np.nanmax(field)), 1.0)
            im = ax.imshow(r, cmap=cmap, vmin=-1.0, vmax=vmax,
                           extent=[-180, 180, -90, 90])
        else:
            vmin, vmax = limits if limits is not None else (float(np.nanmin(field)),
                                                            float(np.nanmax(field)))
            im = ax.imshow(r, cmap=cmap, vmin=vmin, vmax=vmax,
                           extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.72)

    p = outdir / "geography_primitives.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_ocean_geography(world, outdir: Path) -> list[Path]:
    grid = world.grid
    out: list[Path] = []
    basin = world.fields.get("ocean.basin_id")
    depth_province = world.fields.get("ocean.depth_province")
    margin_type = world.fields.get("ocean.margin_type")
    gateway = world.fields.get("ocean.gateway_id")
    elev = world.fields.get("terrain.elevation_m")

    if basin is not None:
        basin = np.asarray(basin, dtype=np.float64)
        if basin.shape == (grid.n,):
            raster = to_raster(grid, basin)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(raster, cmap="tab20", vmin=-1.0,
                           vmax=max(float(np.nanmax(basin)), 1.0),
                           extent=[-180, 180, -90, 90])
            ax.set_title("Ocean basins")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "ocean_basins.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    if depth_province is not None:
        depth_province = np.asarray(depth_province, dtype=np.float64)
        if depth_province.shape == (grid.n,):
            colors = [
                "#efe8d7",  # land
                "#9bd4d8",  # shelf
                "#4aa6b5",  # slope
                "#287c9d",  # rise
                "#173b57",  # abyss
                "#c7b86a",  # ridge
                "#301934",  # trench
                "#d98ca4",  # restricted
            ]
            cmap = ListedColormap(colors)
            norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
            raster = to_raster(grid, depth_province)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
            ax.set_title("Ocean depth provinces")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "ocean_depth_provinces.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    if margin_type is not None:
        margin_type = np.asarray(margin_type, dtype=np.float64)
        if margin_type.shape == (grid.n,):
            colors = [
                "#efe8d7",  # land
                "#7fcdbb",  # passive margin
                "#f03b20",  # active margin
                "#fdbb84",  # ridge
                "#c994c7",  # restricted
                "#225ea8",  # open ocean
            ]
            cmap = ListedColormap(colors)
            norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
            raster = to_raster(grid, margin_type)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
            ax.set_title("Continental and ocean margins")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "continental_margins.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    if elev is not None:
        elev = np.asarray(elev, dtype=np.float64)
        if elev.shape == (grid.n,):
            depth_m = np.where(elev < world.sea_level, world.sea_level - elev, np.nan)
            raster = to_raster(grid, depth_m)
            cmap = plt.get_cmap("cividis_r").copy()
            cmap.set_bad("#efe8d7")
            vmax = max(float(np.nanpercentile(depth_m, 98)) if np.isfinite(depth_m).any()
                       else 1.0, 1000.0)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(np.ma.masked_invalid(raster), cmap=cmap, vmin=0.0, vmax=vmax,
                           extent=[-180, 180, -90, 90])
            ax.set_title("Bathymetry: shelf, slope, abyss")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "bathymetry_shelf_slope_abyss.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    if gateway is not None:
        gateway = np.asarray(gateway, dtype=np.float64)
        if gateway.shape == (grid.n,) and np.any(gateway >= 0.0):
            raster = to_raster(grid, gateway)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(raster, cmap="tab20", vmin=-1.0,
                           vmax=max(float(np.nanmax(gateway)), 1.0),
                           extent=[-180, 180, -90, 90])
            ax.set_title("Ocean gateway components")
            fig.colorbar(im, ax=ax, shrink=0.7)
            p = outdir / "ocean_gateways.png"
            fig.savefig(p, dpi=110, bbox_inches="tight")
            plt.close(fig)
            out.append(p)

    return out


def _render_seam_continuity(world, outdir: Path) -> Path | None:
    from aevum import validation

    detail = validation._seam_continuity_metrics(world)
    metrics = [
        ("land/ocean seam edges", detail["seam_land_ocean_mismatch_fraction"]),
        ("exposed land id seam", detail["seam_exposed_land_component_mismatch_fraction"]),
        ("ocean basin id seam", detail["seam_ocean_basin_mismatch_fraction"]),
        ("render duplicate land", detail["render_duplicate_land_mismatch_fraction"]),
        ("edge-band land", detail["edge_band_land_mismatch_fraction"]),
        ("elev jump ratio", min(detail["seam_to_global_elevation_jump_ratio"] / 3.0, 1.0)),
    ]
    labels = [m[0] for m in metrics]
    values = [float(m[1]) for m in metrics]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    ax = axes[0]
    y = np.arange(len(values))
    ax.barh(y, values, color="#4f83a8")
    ax.set_yticks(y, labels)
    ax.set_xlim(0.0, 1.0)
    ax.set_title("Antimeridian seam diagnostics")
    for yy, val in zip(y, values):
        ax.text(val + 0.01, yy, f"{val:.3f}", va="center", fontsize=8)

    grid = world.grid
    rel = world.get_field("terrain.elevation_m", 0.0) - world.sea_level
    land = rel >= 0.0
    edges = grid.edges
    seam_edges = edges[np.abs(grid.lon[edges[:, 0]] - grid.lon[edges[:, 1]]) > 180.0]
    ax = axes[1]
    ax.set_title("Spherical graph edges crossing seam")
    ax.set_xlim(-181.0, 181.0)
    ax.set_ylim(-90.0, 90.0)
    ax.axvline(-180.0, color="#555555", linewidth=0.8)
    ax.axvline(180.0, color="#555555", linewidth=0.8)
    if seam_edges.size:
        i, j = seam_edges[:, 0], seam_edges[:, 1]
        mismatch = land[i] != land[j]
        lat_mid = 0.5 * (grid.lat[i] + grid.lat[j])
        lon_left = np.where(grid.lon[i] < 0.0, grid.lon[i], grid.lon[j])
        lon_right = np.where(grid.lon[i] >= 0.0, grid.lon[i], grid.lon[j])
        ax.scatter(lon_left[~mismatch], lat_mid[~mismatch], s=8, color="#2c7fb8",
                   alpha=0.75, label="same land/sea")
        ax.scatter(lon_right[~mismatch], lat_mid[~mismatch], s=8, color="#2c7fb8",
                   alpha=0.75)
        ax.scatter(lon_left[mismatch], lat_mid[mismatch], s=12, color="#d7301f",
                   alpha=0.85, label="land/sea change")
        ax.scatter(lon_right[mismatch], lat_mid[mismatch], s=12, color="#d7301f",
                   alpha=0.85)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.legend(loc="lower left", fontsize=7, frameon=False)

    p = outdir / "seam_continuity.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_morphology_diagnostics(world, outdir: Path) -> list[Path]:
    from aevum.diagnostics.morphology import ensure_morphology_fields

    ensure_morphology_fields(world)
    grid = world.grid
    specs = [
        ("tectonics.continental_component_id", "continent_components.png",
         "Continental crust components", "tab20", (-1.0, None)),
        ("tectonics.land_width_steps", "land_width.png",
         "Exposed land width (graph steps)", "viridis", (0.0, None)),
        ("tectonics.continent_width_steps", "continent_width.png",
         "Continental crust width (graph steps)", "viridis", (0.0, None)),
        ("tectonics.ribbon_index", "ribbon_index.png",
         "Exposed land ribbon index", "magma", (0.0, 1.0)),
        ("tectonics.coastline_complexity", "coastline_complexity.png",
         "Exposed-land component complexity", "plasma", (0.0, 1.0)),
        ("tectonics.narrow_neck_index", "narrow_necks.png",
         "Narrow articulation necks", "Reds", (0.0, 1.0)),
    ]
    out: list[Path] = []

    for field_name, filename, title, cmap, limits in specs:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field, dtype=np.float64)
        if field.shape != (grid.n,):
            continue
        raster = to_raster(grid, field)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        if limits[1] is None:
            vmax = float(np.nanpercentile(field[field >= 0.0], 98)) if (field >= 0.0).any() else 1.0
            vmax = max(vmax, 1.0)
        else:
            vmax = float(limits[1])
        im = ax.imshow(raster, cmap=cmap, vmin=float(limits[0]), vmax=vmax,
                       extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.7)
        p = outdir / filename
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        out.append(p)

    contact_specs = [
        ("tectonics.continental_component_id", "Continental components", "tab20", (-1.0, None)),
        ("tectonics.land_width_steps", "Land width", "viridis", (0.0, None)),
        ("tectonics.continent_width_steps", "Continental width", "viridis", (0.0, None)),
        ("tectonics.ribbon_index", "Ribbon index", "magma", (0.0, 1.0)),
        ("tectonics.coastline_complexity", "Coastline complexity", "plasma", (0.0, 1.0)),
        ("tectonics.narrow_neck_index", "Narrow necks", "Reds", (0.0, 1.0)),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(10, 8.2), constrained_layout=True)
    any_panel = False
    for ax, (field_name, title, cmap, limits) in zip(axes.ravel(), contact_specs):
        field = world.fields.get(field_name)
        if field is None or np.asarray(field).shape != (grid.n,):
            ax.axis("off")
            continue
        field = np.asarray(field, dtype=np.float64)
        raster = to_raster(grid, field)
        if limits[1] is None:
            vmax = float(np.nanpercentile(field[field >= 0.0], 98)) if (field >= 0.0).any() else 1.0
            vmax = max(vmax, 1.0)
        else:
            vmax = float(limits[1])
        im = ax.imshow(raster, cmap=cmap, vmin=float(limits[0]), vmax=vmax,
                       extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.72)
        any_panel = True
    p = outdir / "morphology_diagnostics.png"
    if any_panel:
        fig.savefig(p, dpi=110, bbox_inches="tight")
        out.append(p)
    plt.close(fig)
    return out


def _render_tectonic_object_layers(world, outdir: Path) -> list[Path]:
    if "tectonics.continent_width_steps" not in world.fields:
        from aevum.diagnostics.morphology import ensure_morphology_fields
        ensure_morphology_fields(world)
    grid = world.grid
    specs = [
        ("crust.domain", "crust_domain.png", "Crust domain", "tab20", (0.0, 8.0)),
        ("tectonics.continent_id", "tectonic_continent_id.png",
         "Tectonic continent id", "tab20", (-1.0, None)),
        ("tectonics.terrane_id", "terrane_id.png", "Terrane id", "tab20", (-1.0, None)),
        ("tectonics.continent_width_steps", "continent_width.png",
         "Continental crust width (graph steps)", "viridis", (0.0, None)),
        ("tectonics.continental_ribbon_index", "continental_ribbon_index.png",
         "Continental crust ribbon index", "magma", (0.0, 1.0)),
        ("tectonics.continental_narrow_neck_index", "continental_narrow_necks.png",
         "Continental narrow necks", "Reds", (0.0, 1.0)),
    ]
    if not any(name in world.fields for name, *_ in specs):
        return []
    out: list[Path] = []
    for field_name, filename, title, cmap, limits in specs[:3]:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field, dtype=np.float64)
        if field.shape != (grid.n,):
            continue
        raster = to_raster(grid, field)
        if limits[1] is None:
            vmax = float(np.nanmax(field)) if field.size else 1.0
            vmax = max(vmax, 1.0)
        else:
            vmax = float(limits[1])
        fig, ax = plt.subplots(figsize=(9, 4.5))
        im = ax.imshow(raster, cmap=cmap, vmin=float(limits[0]), vmax=vmax,
                       extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.7)
        p = outdir / filename
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        out.append(p)

    fig, axes = plt.subplots(3, 2, figsize=(10, 8.2), constrained_layout=True)
    any_panel = False
    for ax, (field_name, title, cmap, limits) in zip(
        axes.ravel(),
        [(s[0], s[2], s[3], s[4]) for s in specs],
    ):
        field = world.fields.get(field_name)
        if field is None or np.asarray(field).shape != (grid.n,):
            ax.axis("off")
            continue
        field = np.asarray(field, dtype=np.float64)
        raster = to_raster(grid, field)
        if limits[1] is None:
            valid = field[field >= 0.0]
            vmax = float(np.nanmax(valid)) if valid.size else 1.0
            vmax = max(vmax, 1.0)
        else:
            vmax = float(limits[1])
        im = ax.imshow(raster, cmap=cmap, vmin=float(limits[0]), vmax=vmax,
                       extent=[-180, 180, -90, 90])
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.72)
        any_panel = True
    p = outdir / "tectonic_objects.png"
    if any_panel:
        fig.savefig(p, dpi=110, bbox_inches="tight")
        out.append(p)
    plt.close(fig)
    return out


def _render_terrain_province(world, outdir: Path) -> Path | None:
    grid = world.grid
    field = world.fields.get("terrain.province")
    if field is None:
        return None
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (grid.n,):
        return None
    colors = [
        "#173b57",  # ocean
        "#6b8f4e",  # craton/shield
        "#4f9a52",  # continental platform
        "#9fc28b",  # margin/coastal plain
        "#c99454",  # terrane / arc
        "#b35d4d",  # suture/orogen
        "#8f6fb2",  # LIP/plateau
        "#efe8a0",  # highland
        "#78a9c7",  # oceanic island arc
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    raster = to_raster(grid, field)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Terrain provinces")
    fig.colorbar(im, ax=ax, shrink=0.7)
    p = outdir / "terrain_provinces.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_continental_detail(world, outdir: Path) -> Path | None:
    grid = world.grid
    raw_field = world.fields.get("terrain.continental_detail")
    region_field = world.fields.get("terrain.continental_detail_region_code")
    field = region_field if region_field is not None else raw_field
    if field is None:
        return None
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (grid.n,):
        return None
    colors = [
        "#173b57",  # ocean
        "#6f8e59",  # craton shield
        "#4f9a52",  # platform / interior
        "#b7c88a",  # sedimentary basin
        "#87b6a6",  # rift/passive basin
        "#b35d4d",  # orogen
        "#8f6fb2",  # plateau / highland
        "#c99454",  # island arc / microcontinent
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    if region_field is not None and raw_field is not None:
        raw = np.asarray(raw_field, dtype=np.float64)
        if raw.shape == (grid.n,):
            raw_raster = to_raster(grid, raw)
            fig, ax = plt.subplots(figsize=(9, 4.5))
            im = ax.imshow(raw_raster, cmap=cmap, norm=norm,
                           extent=[-180, 180, -90, 90])
            ax.set_title("Continental detail provinces (raw cargo response)")
            fig.colorbar(im, ax=ax, shrink=0.7)
            raw_path = outdir / "continental_detail_raw_provinces.png"
            fig.savefig(raw_path, dpi=110, bbox_inches="tight")
            plt.close(fig)
    raster = to_raster(grid, field)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    title = (
        "Continental detail provinces (regional expression)"
        if region_field is not None else "Continental detail provinces"
    )
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.7)
    p = outdir / "continental_detail_provinces.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_inland_geomorphology_regions(world, outdir: Path) -> Path | None:
    grid = world.grid
    field = world.fields.get("terrain.inland_geomorphology_region_code")
    if field is None:
        return None
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (grid.n,):
        return None
    colors = [
        "#173b57",  # none/ocean
        "#6f8e59",  # shield
        "#5f9f4a",  # platform
        "#b7c88a",  # sag basin
        "#b35d4d",  # old orogen
        "#87b6a6",  # rift
        "#c7aa58",  # platform swell
        "#d98a34",  # escarpment
        "#8f6fb2",  # plateau margin
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    raster = to_raster(grid, field)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Inland geomorphology regions")
    fig.colorbar(im, ax=ax, shrink=0.7)
    p = outdir / "inland_geomorphology_regions.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_internal_geographic_blocks(world, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    specs = (
        ("tectonics.internal_geographic_block_code",
         "internal_geographic_blocks.png",
         "Internal geographic blocks (truth cargo)"),
        ("terrain.internal_geographic_block_region_code",
         "internal_geographic_block_regions.png",
         "Internal geographic blocks (regional expression)"),
    )
    colors = [
        "#173b57",  # none/ocean
        "#6f8e59",  # craton core
        "#a7b76f",  # stable platform
        "#b7c88a",  # intracratonic basin
        "#b35d4d",  # mobile belt
        "#87b6a6",  # rifted margin
        "#c99454",  # accreted terrane
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    for field_name, filename, title in specs:
        field = world.fields.get(field_name)
        if field is None:
            continue
        field = np.asarray(field, dtype=np.float64)
        if field.shape != (world.grid.n,):
            continue
        raster = to_raster(world.grid, field)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
        ax.set_title(title)
        cb = fig.colorbar(im, ax=ax, shrink=0.7)
        cb.set_ticks(list(range(len(colors))))
        cb.set_ticklabels([
            INTERNAL_BLOCK_NAMES.get(i, str(i)) for i in range(len(colors))
        ])
        p = outdir / filename
        fig.savefig(p, dpi=110, bbox_inches="tight")
        plt.close(fig)
        paths.append(p)
    return paths


def _render_wilson_cycle_phase(world, outdir: Path) -> Path | None:
    grid = world.grid
    field = world.fields.get("archive.wilson_cycle_phase")
    if field is None:
        return None
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (grid.n,):
        return None
    colors = [
        "#102a43",  # inactive
        "#8dd3c7",  # rift
        "#4daf4a",  # opening
        "#377eb8",  # mature ocean
        "#ffb000",  # closing
        "#e41a1c",  # collision
        "#7b3294",  # suture
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    raster = to_raster(grid, field)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Wilson-cycle phase")
    fig.colorbar(im, ax=ax, shrink=0.7)
    p = outdir / "wilson_cycle_phase.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def _render_ocean_gateways(world, outdir: Path) -> Path | None:
    grid = world.grid
    gateways = world.objects.get("tectonics.ocean_gateways", [])
    if "ocean.gateway_id" in world.fields:
        return None
    if not gateways:
        return None
    field = np.zeros(grid.n, dtype=np.float64)
    for g in gateways:
        cells = np.asarray(g.get("cells", []), dtype=int)
        cells = cells[(0 <= cells) & (cells < grid.n)]
        if cells.size:
            field[cells] = float(g.get("phase_code", 0.0))
    if not np.any(field > 0.0):
        return None
    colors = [
        "#102a43",
        "#8dd3c7",
        "#4daf4a",
        "#377eb8",
        "#ffb000",
        "#e41a1c",
        "#7b3294",
    ]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), cmap.N)
    raster = to_raster(grid, field)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title("Ocean gateways from Wilson cycles")
    fig.colorbar(im, ax=ax, shrink=0.7)
    p = outdir / "ocean_gateways.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def render_archive_continuity(engine, outdir: Path) -> Path:
    from aevum import validation

    detail = validation.tectonic_diagnostics(engine)["frame_continuity"]
    specs = [
        ("late_max_plate_area_distribution_delta", "plate area"),
        ("late_max_continental_fraction_delta", "continental area"),
        ("late_max_plate_crust_composition_delta", "plate/crust mix"),
        ("late_max_plate_crust_composition_delta_per_100myr", "plate/crust rate"),
        ("late_max_crust_domain_change_fraction", "crust domain"),
        ("late_max_terrain_province_change_fraction", "terrain province"),
        ("late_max_exposed_land_change_fraction", "exposed land"),
        ("late_max_wilson_phase_change_fraction", "Wilson phase"),
    ]
    labels = [label for key, label in specs if key in detail]
    values = [float(detail[key]) for key, _ in specs if key in detail]
    if not values:
        labels = ["frame pairs"]
        values = [float(detail.get("n_frame_pairs", 0))]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    y = np.arange(len(values))
    ax.barh(y, values, color="#4f83a8")
    ax.set_yticks(y, labels)
    ax.set_xlabel("fraction or rate")
    ax.set_title("Archive continuity diagnostics")
    ax.set_xlim(0.0, max(max(values) * 1.15, 0.1))
    for yy, val in zip(y, values):
        ax.text(val, yy, f" {val:.3f}", va="center", fontsize=8)
    p = outdir / "archive_continuity.png"
    fig.savefig(p, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return p


def render_hexmap(cm: CompiledMap, outdir: Path) -> Path:
    def wrapped_columns(arr):
        return np.concatenate([arr[:, -1:], arr, arr[:, :1]], axis=1)

    def wrapped_points(mask):
        rr, cc = np.where(mask)
        rr_all = [rr]
        cc_all = [cc + 1]
        left = cc == 0
        right = cc == cm.hexgrid.width - 1
        if left.any():
            rr_all.append(rr[left])
            cc_all.append(np.full(int(left.sum()), cm.hexgrid.width + 1))
        if right.any():
            rr_all.append(rr[right])
            cc_all.append(np.zeros(int(right.sum()), dtype=int))
        return np.concatenate(rr_all), np.concatenate(cc_all)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.imshow(wrapped_columns(_hex_display_rgb(cm)), aspect="equal")
    ax.axvline(0.5, color="white", linewidth=0.5, alpha=0.7)
    ax.axvline(cm.hexgrid.width + 0.5, color="white", linewidth=0.5, alpha=0.7)

    rr, cc = wrapped_points(cm.river)
    ax.scatter(cc, rr, s=2, c="#1565c0", marker="s")
    # show resources only where they are reachable land/coast tiles
    res_mask = (cm.resources != "") & (cm.terrain >= 1) & (cm.terrain <= 4)
    rr, cc = wrapped_points(res_mask)
    ax.scatter(cc, rr, s=14, c="black", marker="^", edgecolors="white", linewidths=0.3)
    for (r, c) in cm.starts:
        xs = [c + 1]
        if c == 0:
            xs.append(cm.hexgrid.width + 1)
        elif c == cm.hexgrid.width - 1:
            xs.append(0)
        for x in xs:
            ax.scatter([x], [r], s=80, facecolors="none", edgecolors="red", linewidths=2)
            ax.text(x, r, " start", color="red", fontsize=7, va="center")

    handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#173b57",
                   markersize=10, label="ocean"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#5f9fb1",
                   markersize=10, label="coast"),
    ]
    handles.extend(
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=col,
                   markersize=10, label=name)
        for name, col in HEX_LAND_ELEVATION_COLORS[::2]
    )
    handles.append(plt.Line2D([0], [0], marker="s", color="w",
                              markerfacecolor="#eef5f8", markersize=10,
                              label="ice"))
    ax.legend(handles=handles, loc="lower left", fontsize=7, ncol=3)
    ax.set_title("Compiled strategy hex map (rivers, resources ^, starts o)")
    ax.set_xlim(-0.5, cm.hexgrid.width + 1.5)
    ax.set_xticks([]); ax.set_yticks([])
    p = outdir / "hexmap.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    return p


def render_compiler_consistency(cm: CompiledMap, outdir: Path) -> Path | None:
    if cm.source_land_fraction is None or cm.source_shelf_fraction is None:
        return None
    terrain = cm.terrain
    source_land = np.asarray(cm.source_land_fraction, dtype=np.float64)
    source_shelf = np.asarray(cm.source_shelf_fraction, dtype=np.float64)
    compiled_land = (terrain >= 2) & (terrain <= 5)
    broad_land_to_water = (source_land >= 0.70) & ~compiled_land
    broad_ocean_to_land = (source_land <= 0.10) & compiled_land
    shelf_to_deep = (source_shelf >= 0.55) & (source_land < 0.45) & (terrain == 0)
    contradiction = np.zeros(terrain.shape, dtype=np.float64)
    contradiction[broad_land_to_water] = 1.0
    contradiction[broad_ocean_to_land] = 2.0
    contradiction[shelf_to_deep] = 3.0

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), constrained_layout=True)
    specs = [
        (source_land, "Source land fraction", PRECIP_CMAP, (0.0, 1.0)),
        (source_shelf, "Source shelf/shallow fraction", PRECIP_CMAP, (0.0, 1.0)),
        (terrain.astype(float), "Compiled terrain", ListedColormap(TERRAIN_COLORS),
         (-0.5, len(TERRAIN_COLORS) - 0.5)),
        (contradiction, "Compiler contradictions", ListedColormap([
            "#f7f7f7", "#d7301f", "#fc8d59", "#2b8cbe",
        ]), (-0.5, 3.5)),
    ]
    for ax, (arr, title, cmap, limits) in zip(axes.ravel(), specs):
        im = ax.imshow(arr, cmap=cmap, vmin=limits[0], vmax=limits[1], aspect="equal")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.72)
    p = outdir / "compiler_consistency.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return p


def _hex_display_rgb(cm: CompiledMap) -> np.ndarray:
    """Display hexes with ocean/coast classes and a land elevation colour ramp."""
    sea = float(cm.meta.get("sea_level_m", 0.0))
    rel = cm.elevation - sea
    terrain = cm.terrain
    rgb = np.empty(terrain.shape + (3,), dtype=np.float64)

    def color(hex_color: str) -> np.ndarray:
        hex_color = hex_color.lstrip("#")
        return np.array([int(hex_color[i:i + 2], 16) / 255.0
                         for i in (0, 2, 4)], dtype=np.float64)

    rgb[:] = color("#173b57")
    rgb[terrain == 1] = color("#5f9fb1")
    land = (terrain >= 2) & (terrain <= 4)
    thresholds = [0.0, 120.0, 350.0, 700.0, 1100.0, 1700.0, 2500.0, 3600.0]
    colors = [color(c) for _, c in HEX_LAND_ELEVATION_COLORS]
    for k, lo in enumerate(thresholds):
        hi = thresholds[k + 1] if k + 1 < len(thresholds) else np.inf
        mask = land & (rel >= lo) & (rel < hi)
        rgb[mask] = colors[k]
    rgb[terrain == 5] = color("#eef5f8")
    return rgb


def render_timeline(timeline: list[dict], outdir: Path) -> Path:
    types = sorted({e["type"] for e in timeline})
    ymap = {t: i for i, t in enumerate(types)}
    fig, ax = plt.subplots(figsize=(11, 4))
    for e in timeline:
        ax.scatter(e["time_myr"], ymap[e["type"]], s=10 + min(e["magnitude"], 100),
                   alpha=0.6, c="C0")
    ax.set_yticks(range(len(types)))
    ax.set_yticklabels(types, fontsize=8)
    ax.set_xlabel("time since formation (Myr)")
    ax.set_title("Deep-time event timeline")
    ax.grid(True, axis="x", alpha=0.3)
    p = outdir / "timeline.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    return p


def render_history(archive, outdir: Path) -> Path:
    frames = [f for f in archive.frames if "terrain.elevation_m" in f.fields]
    frames = frames[:: max(1, len(frames) // 6)][:6]
    if not frames:
        return outdir / "history.png"
    grid = archive.world.grid
    fig, axes = plt.subplots(2, 3, figsize=(12, 5))
    for ax, fr in zip(axes.ravel(), frames):
        sea = fr.globals.get("ocean.sea_level_m", 0.0)
        r = to_raster(grid, fr.fields["terrain.elevation_m"] - sea, 180, 90)
        ax.imshow(r, cmap=ELEVATION_CMAP, norm=ELEVATION_NORM)
        ax.set_title(f"{fr.time_myr:.0f} Myr", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes.ravel()[len(frames):]:
        ax.axis("off")
    fig.suptitle("Paleogeography through deep time")
    p = outdir / "history.png"
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return p
