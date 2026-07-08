"""Real-Earth climate and biome calibration references.

This module is the first executable piece of the Earth calibration track.  It
keeps source metadata, local cache status, and low-resolution derived products
separate from the generated-world pipeline.  The only direct raster adapter in
this first pass is ETOPO5, because it is already cached locally and can be read
without adding heavy geospatial dependencies.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
import hashlib
import http.client
from io import BytesIO
import json
from pathlib import Path
import re
import tempfile
import time
from typing import Any
from concurrent.futures import ThreadPoolExecutor
import urllib.error
import urllib.parse
import urllib.request
import warnings
import zipfile

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.path import Path as MplPath
from PIL import Image, ImageDraw, ImageFont

from aevum import render
from aevum.core.grid import SphereGrid
from aevum.core.units import CONSTANTS


SCHEMA = "aevum.earth_climate_reference.v1"
SOURCE_MANIFEST_SCHEMA = "aevum.earth_climate_source_manifest.v1"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_ROOT = PROJECT_ROOT / "data" / "reference" / "earth_climate"
ETOPO5_PATH = PROJECT_ROOT / "data" / "reference" / "etopo5" / "ETOPO5.DAT"
ETOPO5_SHAPE = (2160, 4320)
ETOPO2022_DIR = REFERENCE_ROOT / "raw" / "etopo_2022"
ETOPO2022_OPENDAP_ASC = (
    ETOPO2022_DIR / "etopo_2022_v1_60s_bed_0p5deg_opendap.asc"
)
ETOPO2022_OPENDAP_ASC_URL = (
    "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/"
    "60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc.ascii?"
    "lat[0:30:10799],lon[0:30:21599],z[0:30:10799][0:30:21599]"
)
WORLDCLIM_DIR = REFERENCE_ROOT / "raw" / "worldclim_2_1"
WORLDCLIM_TAVG_ZIP = WORLDCLIM_DIR / "wc2.1_10m_tavg.zip"
WORLDCLIM_PREC_ZIP = WORLDCLIM_DIR / "wc2.1_10m_prec.zip"
KOPPEN_ASCII_ZIP = (
    REFERENCE_ROOT
    / "raw"
    / "koppen_geiger"
    / "Koeppen-Geiger-ASCII.zip"
)
GLOH2O_KOPPEN_ZIP = (
    REFERENCE_ROOT
    / "raw"
    / "koppen_geiger"
    / "gloh2o"
    / "koppen_geiger_tif.zip"
)
GLOH2O_KOPPEN_PERIOD = "1991_2020"
GLOH2O_KOPPEN_RESOLUTION = "0p1"
RESOLVE_ECOREGIONS_ZIP = (
    REFERENCE_ROOT
    / "raw"
    / "ecoregions_2017"
    / "Ecoregions2017.zip"
)
NOAA_PSL_NCEP_DIR = REFERENCE_ROOT / "raw" / "noaa_psl_ncep"
NOAA_PSL_UWND_10M = NOAA_PSL_NCEP_DIR / "uwnd.10m.mon.mean.nc"
NOAA_PSL_VWND_10M = NOAA_PSL_NCEP_DIR / "vwnd.10m.mon.mean.nc"
NOAA_PSL_SLP = NOAA_PSL_NCEP_DIR / "slp.mon.mean.nc"
NOAA_OISST_DIR = REFERENCE_ROOT / "raw" / "noaa_oisst_v2"
NOAA_OISST_SST_LTM = NOAA_OISST_DIR / "sst.ltm.1991-2020.nc"
NOAA_OISST_ICEC_LTM = NOAA_OISST_DIR / "icec.ltm.1991-2020.nc"
NOAA_AOML_DRIFTER_DIR = REFERENCE_ROOT / "raw" / "noaa_aoml_drifter_current_v3"
NOAA_AOML_DRIFTER_BASE_URL = (
    "https://tiledimageservices.arcgis.com/P3ePLMYs2RVChkJx/arcgis/rest/services/"
    "annual_drifter_mean_v3/ImageServer"
)
NOAA_AOML_DRIFTER_LEVEL = 2
NOAA_AOML_DRIFTER_EXTENT = (-180.0, -73.0, 180.0, 85.0)
NOAA_AOML_DRIFTER_ORIGIN = (-180.0, 85.0)
NOAA_AOML_DRIFTER_RESOLUTIONS = {0: 1.0, 1: 0.5, 2: 0.25}
OSCAR_ARCGIS_DIR = REFERENCE_ROOT / "raw" / "oscar" / "monthly_mean_2001_2020_arcgis"
OSCAR_ARCGIS_ITEM_URL = (
    "https://www.arcgis.com/home/item.html?id=b02f417ebbed4dc69edefd848dc69715"
)
OSCAR_ARCGIS_SERVICE_URL = (
    "https://tiledimageservices.arcgis.com/P3ePLMYs2RVChkJx/arcgis/rest/services/"
    "OSCAR_2001_to_2020_20yr_Monthly_Mean/ImageServer"
)
OSCAR_ARCGIS_MULTIDIMENSIONAL_INFO = OSCAR_ARCGIS_DIR / "multidimensional_info.json"
OSCAR_ARCGIS_LEVEL = 1
OSCAR_ARCGIS_EXTENT = (
    -20037507.0672,
    -30240971.45838615,
    20051950.036327545,
    30240971.45838615,
)
OSCAR_ARCGIS_ORIGIN = (-20037507.0672, 30240971.4583861)
OSCAR_ARCGIS_RESOLUTIONS = {
    0: 222564.647347828,
    1: 111282.323673914,
    2: 55641.161836957,
    3: 27820.5809184785,
}
OSCAR_WEBMERCATOR_RADIUS_M = 6378137.0
OSCAR_WEBMERCATOR_MAX_LAT_DEG = 85.05112878
ESA_CCI_LAND_COVER_DIR = (
    REFERENCE_ROOT / "raw" / "esa_cci_land_cover" / "planetary_computer_preview_2020"
)
ESA_CCI_LAND_COVER_COLLECTION_URL = (
    "https://planetarycomputer.microsoft.com/api/stac/v1/collections/esa-cci-lc"
)
ESA_CCI_LAND_COVER_ITEMS_URL = f"{ESA_CCI_LAND_COVER_COLLECTION_URL}/items?limit=100"
ESA_CCI_LAND_COVER_YEAR = 2020
ESA_CCI_LAND_COVER_PREVIEW_SIZE = 1024
ESA_CCI_LAND_COVER_TILE_DEG = 45.0
ARCGIS_TILE_SIZE = 256

SEASON_LABELS = ("DJF", "MAM", "JJA", "SON")
SEASON_MONTHS = ((11, 0, 1), (2, 3, 4), (5, 6, 7), (8, 9, 10))
KOPPEN_CLASSES = (
    "Af", "Am", "As", "Aw",
    "BSh", "BSk", "BWh", "BWk",
    "Cfa", "Cfb", "Cfc", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc",
    "Dfa", "Dfb", "Dfc", "Dfd", "Dsa", "Dsb", "Dsc", "Dsd",
    "Dwa", "Dwb", "Dwc", "Dwd",
    "EF", "ET",
)
KOPPEN_TO_CODE = {name: idx + 1 for idx, name in enumerate(KOPPEN_CLASSES)}
KOPPEN_CODE_TO_CLASS = {idx + 1: name for idx, name in enumerate(KOPPEN_CLASSES)}
GLOH2O_KOPPEN_CODE_TO_CLASS = {
    1: "Af", 2: "Am", 3: "Aw",
    4: "BWh", 5: "BWk", 6: "BSh", 7: "BSk",
    8: "Csa", 9: "Csb", 10: "Csc",
    11: "Cwa", 12: "Cwb", 13: "Cwc",
    14: "Cfa", 15: "Cfb", 16: "Cfc",
    17: "Dsa", 18: "Dsb", 19: "Dsc", 20: "Dsd",
    21: "Dwa", 22: "Dwb", 23: "Dwc", 24: "Dwd",
    25: "Dfa", 26: "Dfb", 27: "Dfc", 28: "Dfd",
    29: "ET", 30: "EF",
}
KOPPEN_MAJOR_COLORS = [
    "#173b57",  # ocean / no-data
    "#176c43",  # A tropical
    "#d9c98f",  # B dry
    "#6ba35d",  # C temperate
    "#9fb0cf",  # D continental
    "#eef5f8",  # E polar
]
KOPPEN_MAJOR_CMAP = ListedColormap(KOPPEN_MAJOR_COLORS)
KOPPEN_MAJOR_NORM = BoundaryNorm(np.arange(-0.5, 6.5, 1.0), KOPPEN_MAJOR_CMAP.N)
GLOH2O_KOPPEN_COLORS = [
    "#173b57",
    "#0000ff", "#0078ff", "#46aafa",
    "#ff0000", "#ff9696", "#f5a500", "#ffdc64",
    "#ffff00", "#c8c800", "#969600",
    "#96ff96", "#64c864", "#329632",
    "#c8ff50", "#64ff50", "#32c800",
    "#ff00ff", "#c800c8", "#963296", "#966496",
    "#aaafff", "#5a78dc", "#4b50b4", "#320087",
    "#00ffff", "#37c8ff", "#007d7d", "#00465f",
    "#b2b2b2", "#686868",
]
GLOH2O_KOPPEN_CMAP = ListedColormap(GLOH2O_KOPPEN_COLORS)
GLOH2O_KOPPEN_NORM = BoundaryNorm(
    np.arange(-0.5, len(GLOH2O_KOPPEN_COLORS) + 0.5, 1.0),
    GLOH2O_KOPPEN_CMAP.N,
)
RESOLVE_BIOME_NAMES = {
    0: "ocean / no-data",
    1: "Tropical & Subtropical Moist Broadleaf Forests",
    2: "Tropical & Subtropical Dry Broadleaf Forests",
    3: "Tropical & Subtropical Coniferous Forests",
    4: "Temperate Broadleaf & Mixed Forests",
    5: "Temperate Conifer Forests",
    6: "Boreal Forests/Taiga",
    7: "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    8: "Temperate Grasslands, Savannas & Shrublands",
    9: "Flooded Grasslands & Savannas",
    10: "Montane Grasslands & Shrublands",
    11: "Tundra",
    12: "Mediterranean Forests, Woodlands & Scrub",
    13: "Deserts & Xeric Shrublands",
    14: "Mangroves",
}
RESOLVE_BIOME_COLORS = [
    "#173b57",
    "#0b6b3a", "#9bc35b", "#238b45", "#63a35c", "#2f7f5f",
    "#7da0a6", "#c4b85f", "#d6c27b", "#5fa4a3", "#8e9e63",
    "#d8e4e8", "#c77948", "#d9c98f", "#1b7f68",
]
RESOLVE_BIOME_CMAP = ListedColormap(RESOLVE_BIOME_COLORS)
RESOLVE_BIOME_NORM = BoundaryNorm(
    np.arange(-0.5, len(RESOLVE_BIOME_COLORS) + 0.5, 1.0),
    RESOLVE_BIOME_CMAP.N,
)
LAND_COVER_BROAD_NAMES = {
    0: "no-data",
    1: "water",
    2: "cropland",
    3: "forest",
    4: "grass-shrub",
    5: "wetland",
    6: "urban",
    7: "bare-sparse",
    8: "snow-ice",
}
LAND_COVER_BROAD_COLORS = [
    "#ffffff",
    "#0046c8",
    "#fff064",
    "#006400",
    "#be9600",
    "#00785a",
    "#c31400",
    "#fff5d7",
    "#f0f0f0",
]
LAND_COVER_BROAD_CMAP = ListedColormap(LAND_COVER_BROAD_COLORS)
LAND_COVER_BROAD_NORM = BoundaryNorm(
    np.arange(-0.5, len(LAND_COVER_BROAD_COLORS) + 0.5, 1.0),
    LAND_COVER_BROAD_CMAP.N,
)


@dataclass(frozen=True)
class EarthClimateReferenceConfig:
    cells: tuple[int, ...] = (8000,)
    width: int = 720
    height: int = 360
    render_assets: bool = True
    manifest_out: Path | None = None
    include_worldclim: bool = True
    include_koppen: bool = True
    include_noaa_psl: bool = True
    include_sst: bool = True
    include_ocean_currents: bool = True
    include_etopo2022: bool = True
    include_land_cover: bool = True
    download_oscar: bool = False
    download_land_cover: bool = False
    climatology_start_year: int = 1991
    climatology_end_year: int = 2020


REFERENCE_SOURCES: tuple[dict[str, Any], ...] = (
    {
        "source_id": "NOAA_ETOPO5_LOCAL",
        "category": "topography_bathymetry",
        "title": "NOAA ETOPO5 5-minute global relief grid",
        "source_url": "https://www.ngdc.noaa.gov/mgg/global/etopo5.HTML",
        "variables": ("elevation_m", "bathymetry_m"),
        "temporal_baseline": "static modern Earth relief",
        "native_resolution": "5 arc-minute latitude/longitude grid",
        "local_cache_path": "data/reference/etopo5/ETOPO5.DAT",
        "requires_account": False,
        "priority": "available_fallback",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.elevation_m",
            "earth.land_mask",
            "earth.hypsometry_metrics",
        ),
        "notes": (
            "Already cached in this workspace.  Used as the first executable "
            "same-grid Earth topography/bathymetry reference."
        ),
    },
    {
        "source_id": "NOAA_ETOPO_2022",
        "category": "topography_bathymetry",
        "title": "NOAA NCEI ETOPO 2022 Global Relief Model",
        "source_url": "https://www.ncei.noaa.gov/products/etopo-global-relief-model",
        "download_urls": (
            "https://www.ngdc.noaa.gov/thredds/catalog/global/ETOPO2022/60s/60s_bed_elev_netcdf/catalog.html",
            ETOPO2022_OPENDAP_ASC_URL,
        ),
        "variables": ("elevation_m", "bathymetry_m"),
        "temporal_baseline": "static modern Earth relief",
        "native_resolution": "15 arc-second source; cached 60 arc-second bedrock relief sampled every 30 cells to 0.5 degree",
        "local_cache_path": "data/reference/earth_climate/raw/etopo_2022/",
        "requires_account": False,
        "priority": "primary_topography_crosscheck",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.etopo2022_elevation_m",
            "earth.etopo2022_land_mask",
            "earth.etopo2022_crosscheck_metrics",
        ),
        "notes": (
            "Used as a same-grid cross-check against the legacy ETOPO5 baseline. "
            "The default pipeline keeps ETOPO5 as earth.elevation_m for backward "
            "compatibility and stores ETOPO 2022 under explicit etopo2022 fields."
        ),
    },
    {
        "source_id": "WORLDCLIM_2_1",
        "category": "temperature_precipitation",
        "title": "WorldClim v2.1 historical monthly climate normals",
        "source_url": "https://www.worldclim.org/data/worldclim21.html",
        "variables": ("monthly_temperature_C", "monthly_precipitation_mm"),
        "temporal_baseline": "1970-2000 normals",
        "native_resolution": "30 arc-second to 10 arc-minute global rasters",
        "local_cache_path": "data/reference/earth_climate/raw/worldclim_2_1/",
        "download_urls": (
            "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_tavg.zip",
            "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_prec.zip",
        ),
        "requires_account": False,
        "priority": "primary_land_climatology",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.monthly_temperature_C",
            "earth.monthly_precip_mm",
            "earth.seasonal_temperature_C",
            "earth.seasonal_precip_mm_yr_equiv",
        ),
        "notes": "Land climatology anchor for seasonal temperature and precipitation.",
    },
    {
        "source_id": "KOTTEK_RUBEL_KOPPEN_2006_ASCII",
        "category": "climate_class",
        "title": "World map of the Koppen-Geiger climate classification updated",
        "source_url": "https://koeppen-geiger.vu-wien.ac.at/present.htm",
        "download_urls": (
            "https://koeppen-geiger.vu-wien.ac.at/data/Koeppen-Geiger-ASCII.zip",
        ),
        "variables": ("koppen_geiger_class",),
        "temporal_baseline": "1951-2000",
        "native_resolution": "0.5 degree ASCII grid",
        "local_cache_path": "data/reference/earth_climate/raw/koppen_geiger/Koeppen-Geiger-ASCII.zip",
        "requires_account": False,
        "priority": "available_climate_class_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.koppen_class",
            "earth.koppen_major_class",
            "earth.biome_class_proxy",
        ),
        "notes": (
            "Used as a lightweight executable climate-class reference.  GloH2O "
            "V3 remains the planned high-resolution follow-up."
        ),
    },
    {
        "source_id": "NOAA_PSL_NCEP_NCAR_REANALYSIS_1",
        "category": "wind_pressure_reanalysis",
        "title": "NOAA PSL NCEP/NCAR Reanalysis 1",
        "source_url": "https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html",
        "download_urls": (
            "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/uwnd.10m.mon.mean.nc",
            "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/vwnd.10m.mon.mean.nc",
            "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.mon.mean.nc",
        ),
        "variables": ("u10", "v10", "surface_pressure", "air_temperature"),
        "temporal_baseline": "1948-present monthly means",
        "native_resolution": "2.5 degree global grids",
        "local_cache_path": "data/reference/earth_climate/raw/noaa_psl_ncep/",
        "requires_account": False,
        "priority": "fallback_global_reanalysis",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.wind_u10_v10",
            "earth.seasonal_pressure_proxy",
        ),
        "notes": "No-account fallback for wind and pressure calibration.",
    },
    {
        "source_id": "ERA5_MONTHLY_SINGLE_LEVELS",
        "category": "wind_pressure_reanalysis",
        "title": "ERA5 monthly averaged single levels",
        "source_url": "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means",
        "variables": ("u10", "v10", "msl", "t2m", "tp"),
        "temporal_baseline": "1940-present monthly means",
        "native_resolution": "global reanalysis grid, user-selected download resolution",
        "local_cache_path": "data/reference/earth_climate/raw/era5/",
        "requires_account": False,
        "priority": "optional_high_quality_reanalysis",
        "parser_status": "planned_optional",
        "derived_targets": (
            "earth.wind_u10_v10",
            "earth.monthly_temperature_C",
            "earth.monthly_precip_mm",
        ),
        "notes": "High-quality reference, but should not block no-account calibration.",
    },
    {
        "source_id": "NOAA_OISST_V2_LTM_1991_2020",
        "category": "sea_surface_temperature",
        "title": "NOAA Optimum Interpolation SST v2 monthly long-term means",
        "source_url": "https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html",
        "download_urls": (
            "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2/sst.ltm.1991-2020.nc",
            "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2/icec.ltm.1991-2020.nc",
        ),
        "variables": ("monthly_sst_C", "monthly_sea_ice_concentration_pct"),
        "temporal_baseline": "1991-2020 monthly long-term mean",
        "native_resolution": "1 degree latitude/longitude grid",
        "local_cache_path": "data/reference/earth_climate/raw/noaa_oisst_v2/",
        "requires_account": False,
        "priority": "primary_sst_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.monthly_sst_C",
            "earth.seasonal_sst_C",
            "earth.annual_sst_C",
            "earth.monthly_sea_ice_concentration_pct",
        ),
        "notes": (
            "Small no-account SST/sea-ice climatology used to calibrate ocean "
            "thermal gradients before full coupled SST dynamics."
        ),
    },
    {
        "source_id": "NOAA_AOML_DRIFTER_CURRENT_CLIMATOLOGY_V3",
        "category": "surface_ocean_currents",
        "title": "NOAA/AOML annual near-surface drifter velocity climatology v3",
        "source_url": "https://www.aoml.noaa.gov/phod/gdp/mean_velocity.php",
        "download_urls": (
            "https://www.arcgis.com/home/item.html?id=3f453a562771441f9d42a2f03c9b6111",
        ),
        "variables": ("annual_current_speed_m_s", "annual_current_direction_deg"),
        "temporal_baseline": "2005-2023 annual climatology",
        "native_resolution": "0.25 degree ArcGIS LERC raster tiles",
        "local_cache_path": "data/reference/earth_climate/raw/noaa_aoml_drifter_current_v3/",
        "requires_account": False,
        "priority": "primary_no_account_current_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.surface_current_u_v",
            "earth.annual_surface_current_speed_m_s",
            "earth.annual_surface_current_direction_deg",
        ),
        "notes": (
            "Public numerical LERC tile service.  Provides annual observed "
            "near-surface current magnitude/direction; OSCAR monthly currents "
            "remain the higher-cost optional follow-up."
        ),
    },
    {
        "source_id": "NASA_JPL_OSCAR",
        "category": "surface_ocean_currents",
        "title": "NASA/JPL PO.DAAC OSCAR surface currents monthly means",
        "source_url": "https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0",
        "download_urls": (
            OSCAR_ARCGIS_ITEM_URL,
            OSCAR_ARCGIS_SERVICE_URL,
            "https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/oscar/open/L4/oscar_v2.0/docs/oscarv2guide.pdf",
        ),
        "variables": (
            "monthly_surface_current_u",
            "monthly_surface_current_v",
            "seasonal_surface_current_u",
            "seasonal_surface_current_v",
            "annual_surface_current_u",
            "annual_surface_current_v",
        ),
        "temporal_baseline": "2001-2020 monthly means sampled as a monthly climatology",
        "native_resolution": "ArcGIS LERC2D Vector-UV cache; source OSCAR v2.0 0.25 degree, default sampled at level 1 for global calibration",
        "local_cache_path": "data/reference/earth_climate/raw/oscar/monthly_mean_2001_2020_arcgis/",
        "requires_account": False,
        "priority": "primary_current_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.monthly_surface_current_u_v",
            "earth.seasonal_surface_current_u_v",
            "earth.annual_surface_current_u_v",
            "earth.monthly_surface_current_speed_m_s",
            "earth.seasonal_surface_current_speed_m_s",
            "earth.annual_surface_current_speed_m_s",
        ),
        "notes": (
            "Public ArcGIS tiled imagery layer derived from OSCAR and credited to "
            "NASA PO.DAAC.  Cache construction is explicit because a full 2001-2020 "
            "monthly climatology requires many LERC tile reads."
        ),
    },
    {
        "source_id": "GLOH2O_KOPPEN_GEIGER",
        "category": "climate_class",
        "title": "Beck et al. / GloH2O Koppen-Geiger climate maps",
        "source_url": "https://www.gloh2o.org/koppen/",
        "download_urls": (
            "https://ndownloader.figshare.com/files/42602809",
        ),
        "variables": ("koppen_geiger_class",),
        "temporal_baseline": "1991-2020 historical classification",
        "native_resolution": "0.01, 0.1, 0.5, and 1.0 degree GeoTIFF rasters",
        "local_cache_path": "data/reference/earth_climate/raw/koppen_geiger/",
        "requires_account": False,
        "priority": "primary_climate_class_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.koppen_class",
            "earth.koppen_major_class",
            "earth.biome_class_proxy",
        ),
        "notes": (
            "Primary executable Koppen source.  The adapter samples the "
            "1991-2020 0.1 degree GeoTIFF by default and falls back to the "
            "older ASCII product only when this cache is unavailable."
        ),
    },
    {
        "source_id": "ESA_CCI_LAND_COVER",
        "category": "biome_land_cover",
        "title": "ESA CCI / C3S global land cover via Planetary Computer",
        "source_url": "https://planetarycomputer.microsoft.com/dataset/esa-cci-lc",
        "download_urls": (
            ESA_CCI_LAND_COVER_COLLECTION_URL,
            ESA_CCI_LAND_COVER_ITEMS_URL,
            "https://doi.org/10.24381/cds.006f2c9a",
        ),
        "variables": ("lccs_class", "broad_land_cover_class"),
        "temporal_baseline": "2020 annual ESA CCI/C3S land cover",
        "native_resolution": (
            "source 300 m COG tiles; cached 1024x1024 categorical rendered "
            "preview per 45 degree tile for coarse calibration"
        ),
        "local_cache_path": (
            "data/reference/earth_climate/raw/esa_cci_land_cover/"
            "planetary_computer_preview_2020/"
        ),
        "requires_account": False,
        "priority": "primary_land_cover_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.esa_cci_lccs_class",
            "earth.esa_cci_land_cover_broad_class",
            "earth.land_cover_crosscheck_metrics",
        ),
        "notes": (
            "Uses Planetary Computer STAC metadata and rendered categorical "
            "previews to avoid downloading full 300 m global COGs.  Suitable "
            "for coarse biome/land-cover envelope checks, not fine land-use "
            "analysis."
        ),
    },
    {
        "source_id": "MODIS_MCD12C1",
        "category": "biome_land_cover",
        "title": "MODIS MCD12C1 global land-cover type",
        "source_url": "https://modis.gsfc.nasa.gov/data/dataprod/mod12.php",
        "variables": ("land_cover_class",),
        "temporal_baseline": "satellite-era annual land cover",
        "native_resolution": "0.05 degree global product",
        "local_cache_path": "data/reference/earth_climate/raw/modis_mcd12c1/",
        "requires_account": False,
        "priority": "secondary_land_cover_reference",
        "parser_status": "planned",
        "derived_targets": ("earth.biome_class",),
        "notes": "Alternative or cross-check for biome envelope calibration.",
    },
    {
        "source_id": "RESOLVE_ECOREGIONS_2017",
        "category": "biome_ecoregion",
        "title": "RESOLVE Ecoregions 2017",
        "source_url": "https://ecoregions.appspot.com/",
        "download_urls": (
            "https://storage.googleapis.com/teow2016/Ecoregions2017.zip",
        ),
        "variables": ("biome_id", "ecoregion_id"),
        "temporal_baseline": "modern terrestrial ecoregions",
        "native_resolution": "global polygon dataset",
        "local_cache_path": "data/reference/earth_climate/raw/ecoregions_2017/",
        "requires_account": False,
        "priority": "secondary_biome_reference",
        "parser_status": "implemented",
        "derived_targets": (
            "earth.resolve_biome_class",
            "earth.resolve_ecoregion_id",
        ),
        "notes": (
            "Terrestrial 14-biome and 846-ecoregion reference.  This is a "
            "calibration target, not a direct climate-driver field."
        ),
    },
)


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _local_path(raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_source_manifest() -> dict[str, Any]:
    """Return source metadata plus current local cache status."""
    entries: list[dict[str, Any]] = []
    for source in REFERENCE_SOURCES:
        entry = dict(source)
        local_path = _local_path(str(source["local_cache_path"]))
        entry["local_cache_abs_path"] = str(local_path)
        entry["local_cache_exists"] = bool(local_path.exists())
        if local_path.is_file():
            entry["local_cache_bytes"] = int(local_path.stat().st_size)
            entry["checksum_sha256"] = _sha256(local_path)
            entry["checksum_status"] = "computed"
        elif local_path.is_dir():
            files = [p for p in local_path.rglob("*") if p.is_file()]
            entry["local_cache_file_count"] = len(files)
            entry["local_cache_bytes"] = int(sum(p.stat().st_size for p in files))
            entry["local_cache_files"] = [
                {
                    "path": str(p.relative_to(PROJECT_ROOT)),
                    "bytes": int(p.stat().st_size),
                    "sha256": _sha256(p),
                }
                for p in sorted(files)[:32]
            ]
            entry["checksum_sha256"] = None
            entry["checksum_status"] = (
                "file_checksums_computed"
                if len(files) <= 32 else "partial_file_checksums_computed"
            )
        else:
            entry["local_cache_bytes"] = 0
            entry["checksum_sha256"] = None
            entry["checksum_status"] = "pending_until_acquisition"
        entry["acquisition_status"] = (
            "available_local_cache" if entry["local_cache_exists"]
            else "not_downloaded"
        )
        entries.append(entry)

    return {
        "schema": SOURCE_MANIFEST_SCHEMA,
        "manifest_version": "2026-07-05.r0",
        "raw_data_policy": (
            "Keep large raw Earth datasets out of git.  Store public source "
            "metadata, acquisition scripts, checksums, and derived low-resolution "
            "fixtures that can be regenerated."
        ),
        "source_count": len(entries),
        "available_local_source_count": int(sum(e["local_cache_exists"] for e in entries)),
        "implemented_parser_count": int(sum(e["parser_status"] == "implemented" for e in entries)),
        "entries": entries,
    }


def write_source_manifest(path: Path = REFERENCE_ROOT / "source_manifest.json") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest = build_source_manifest()
    path.write_text(json.dumps(manifest, indent=2, default=_json_default))
    return path


def load_etopo5(path: Path = ETOPO5_PATH) -> np.ndarray:
    """Load the locally cached ETOPO5 big-endian int16 grid in meters."""
    if not path.exists():
        raise FileNotFoundError(path)
    raw = np.fromfile(path, dtype=">i2")
    expected = ETOPO5_SHAPE[0] * ETOPO5_SHAPE[1]
    if raw.size != expected:
        raise ValueError(f"expected {expected} ETOPO5 cells, found {raw.size}")
    return raw.reshape(ETOPO5_SHAPE).astype(np.float32)


def sample_latlon_raster_to_grid(
    raster: np.ndarray,
    grid: SphereGrid,
    *,
    north_to_south: bool = True,
    lon_origin_deg: float = 0.0,
) -> np.ndarray:
    """Nearest-neighbour sample a regular global lat/lon raster to a sphere grid."""
    raster = np.asarray(raster)
    if raster.ndim != 2:
        raise ValueError("raster must be a 2-D lat/lon array")
    height, width = raster.shape
    lat_step = 180.0 / float(height)
    lon_step = 360.0 / float(width)
    if north_to_south:
        lat_idx = np.rint((90.0 - lat_step / 2.0 - grid.lat) / lat_step)
    else:
        lat_idx = np.rint((grid.lat + 90.0 - lat_step / 2.0) / lat_step)
    lon = (grid.lon - float(lon_origin_deg) + 360.0) % 360.0
    lon_idx = np.rint((lon - lon_step / 2.0) / lon_step)
    lat_idx = np.clip(lat_idx.astype(np.int64), 0, height - 1)
    lon_idx = lon_idx.astype(np.int64) % width
    return raster[lat_idx, lon_idx]


def sample_etopo5_to_grid(grid: SphereGrid, path: Path = ETOPO5_PATH) -> np.ndarray:
    return sample_latlon_raster_to_grid(load_etopo5(path), grid).astype(np.float32)


def etopo2022_available(path: Path = ETOPO2022_OPENDAP_ASC) -> bool:
    return path.exists() and path.stat().st_size > 0


def _parse_opendap_vector(
    text: str,
    name: str,
    *,
    expected: int | None = None,
) -> np.ndarray:
    match = re.search(rf"(?m)^{re.escape(name)}\[(\d+)\]\s*$", text)
    if match is None:
        raise ValueError(f"missing OPeNDAP vector {name!r}")
    count = int(match.group(1))
    start = match.end()
    end = text.find("\n\n", start)
    if end < 0:
        end = len(text)
    values = np.fromstring(text[start:end].replace(",", " "), sep=" ")
    target = count if expected is None else int(expected)
    if values.size != target or count != target:
        raise ValueError(
            f"expected {target} values for {name}, found header {count} and data {values.size}"
        )
    return values.astype(np.float32)


def load_etopo2022_opendap_ascii(
    path: Path = ETOPO2022_OPENDAP_ASC,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load a lightweight OPeNDAP ASCII cache of ETOPO 2022 bedrock relief.

    The cache is produced from NOAA THREDDS with a stride query over the 60
    arc-second global bedrock product.  It is intentionally a small derived
    reference product, not the full 491 MB NetCDF.  Returned latitude is
    ascending south-to-north, matching the OPeNDAP response.
    """
    if not path.exists():
        raise FileNotFoundError(path)
    text = path.read_text(errors="ignore")
    lat = _parse_opendap_vector(text, "lat")
    lon = _parse_opendap_vector(text, "lon")

    marker = re.search(r"(?m)^z\.z\[(\d+)\]\[(\d+)\]\s*$", text)
    if marker is None:
        marker = re.search(r"(?m)^z\[(\d+)\]\[(\d+)\]\s*$", text)
    if marker is None:
        raise ValueError(f"missing ETOPO 2022 z matrix in {path}")
    n_lat = int(marker.group(1))
    n_lon = int(marker.group(2))
    rows: list[np.ndarray] = []
    for line in text[marker.end():].splitlines():
        line = line.strip()
        if not line:
            if rows:
                break
            continue
        if not line.startswith("["):
            if rows:
                break
            continue
        try:
            _, payload = line.split(",", 1)
        except ValueError as exc:
            raise ValueError(f"malformed ETOPO 2022 matrix row: {line[:80]!r}") from exc
        values = np.fromstring(payload.replace(",", " "), sep=" ", dtype=np.float32)
        if values.size != n_lon:
            raise ValueError(
                f"expected {n_lon} ETOPO 2022 columns, found {values.size}"
            )
        rows.append(values)
        if len(rows) == n_lat:
            break
    if len(rows) != n_lat:
        raise ValueError(f"expected {n_lat} ETOPO 2022 rows, found {len(rows)}")
    if lat.size != n_lat or lon.size != n_lon:
        raise ValueError(
            f"ETOPO 2022 coordinate shape mismatch: lat {lat.size}, lon {lon.size}, z {(n_lat, n_lon)}"
        )
    raster = np.vstack(rows).astype(np.float32)
    raster = np.where(np.isfinite(raster), raster, np.nan).astype(np.float32)
    return raster, lat, lon


def sample_etopo2022_to_grid(
    grid: SphereGrid,
    path: Path = ETOPO2022_OPENDAP_ASC,
) -> np.ndarray:
    raster, lat, lon = load_etopo2022_opendap_ascii(path)
    return sample_latlon_coordinate_raster_to_grid(raster, lat, lon, grid).astype(np.float32)


def etopo2022_crosscheck_metrics(
    grid: SphereGrid,
    etopo5_elevation_m: np.ndarray,
    etopo2022_elevation_m: np.ndarray,
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    e5 = np.asarray(etopo5_elevation_m, dtype=np.float64)
    e22 = np.asarray(etopo2022_elevation_m, dtype=np.float64)
    valid = np.isfinite(e5) & np.isfinite(e22)
    if not valid.any():
        return {
            "valid_area_fraction": 0.0,
            "mean_abs_elevation_delta_m": np.nan,
            "p95_abs_elevation_delta_m": np.nan,
            "land_fraction_delta": np.nan,
            "land_elevation_mean_delta_m": np.nan,
            "ocean_depth_mean_delta_m": np.nan,
            "land_ocean_class_mismatch_area_fraction": np.nan,
        }
    diff = e22[valid] - e5[valid]
    weights = area[valid]
    m5 = earth_elevation_metrics(grid, e5)
    m22 = earth_elevation_metrics(grid, e22)
    class_mismatch = valid & ((e5 >= 0.0) != (e22 >= 0.0))
    return {
        "valid_area_fraction": float(area[valid].sum() / total_area),
        "mean_elevation_delta_m": float(np.average(diff, weights=weights)),
        "mean_abs_elevation_delta_m": float(np.average(np.abs(diff), weights=weights)),
        "p50_abs_elevation_delta_m": float(np.percentile(np.abs(diff), 50)),
        "p95_abs_elevation_delta_m": float(np.percentile(np.abs(diff), 95)),
        "land_fraction_delta": float(m22["land_fraction"] - m5["land_fraction"]),
        "land_elevation_mean_delta_m": float(
            m22["land_elevation_mean_m"] - m5["land_elevation_mean_m"]
        ),
        "ocean_depth_mean_delta_m": float(
            m22["ocean_depth_mean_m"] - m5["ocean_depth_mean_m"]
        ),
        "shelf_fraction_of_ocean_delta": float(
            m22["shelf_fraction_of_ocean"] - m5["shelf_fraction_of_ocean"]
        ),
        "abyss_fraction_of_ocean_delta": float(
            m22["abyss_fraction_of_ocean"] - m5["abyss_fraction_of_ocean"]
        ),
        "land_ocean_class_mismatch_area_fraction": float(
            area[class_mismatch].sum() / total_area
        ),
    }


def _read_tif_from_zip(path: Path, member: str) -> np.ndarray:
    with zipfile.ZipFile(path) as zf:
        data = zf.read(member)
    with Image.open(BytesIO(data)) as image:
        return np.asarray(image).astype(np.float32)


def worldclim_available() -> bool:
    return WORLDCLIM_TAVG_ZIP.exists() and WORLDCLIM_PREC_ZIP.exists()


def load_worldclim_monthly(
    *,
    tavg_zip: Path = WORLDCLIM_TAVG_ZIP,
    prec_zip: Path = WORLDCLIM_PREC_ZIP,
) -> tuple[np.ndarray, np.ndarray]:
    """Load WorldClim v2.1 10-minute monthly temperature and precipitation.

    Returns:
        `(monthly_temperature_C, monthly_precip_mm)`, both shaped
        `(12, 1080, 2160)` with no-data cells set to `NaN`.
    """
    if not tavg_zip.exists() or not prec_zip.exists():
        raise FileNotFoundError(
            f"WorldClim zips missing: {tavg_zip} / {prec_zip}"
        )
    tavg: list[np.ndarray] = []
    prec: list[np.ndarray] = []
    with zipfile.ZipFile(tavg_zip) as zf:
        names = set(zf.namelist())
    for month in range(1, 13):
        member = f"wc2.1_10m_tavg_{month:02d}.tif"
        if member not in names:
            raise ValueError(f"missing {member} in {tavg_zip}")
        arr = _read_tif_from_zip(tavg_zip, member)
        arr = np.where(arr < -1.0e20, np.nan, arr)
        tavg.append(arr.astype(np.float32))

    with zipfile.ZipFile(prec_zip) as zf:
        names = set(zf.namelist())
    for month in range(1, 13):
        member = f"wc2.1_10m_prec_{month:02d}.tif"
        if member not in names:
            raise ValueError(f"missing {member} in {prec_zip}")
        arr = _read_tif_from_zip(prec_zip, member)
        arr = np.where(arr <= -32000.0, np.nan, arr)
        prec.append(arr.astype(np.float32))
    return np.asarray(tavg), np.asarray(prec)


def sample_worldclim_to_grid(
    grid: SphereGrid,
) -> dict[str, np.ndarray]:
    monthly_temp_raster, monthly_precip_raster = load_worldclim_monthly()
    monthly_temp = np.vstack([
        sample_latlon_raster_to_grid(
            monthly_temp_raster[i],
            grid,
            north_to_south=True,
            lon_origin_deg=-180.0,
        )[None, :]
        for i in range(12)
    ]).astype(np.float32)
    monthly_precip = np.vstack([
        sample_latlon_raster_to_grid(
            monthly_precip_raster[i],
            grid,
            north_to_south=True,
            lon_origin_deg=-180.0,
        )[None, :]
        for i in range(12)
    ]).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        seasonal_temp = np.vstack([
            np.nanmean(monthly_temp[list(months)], axis=0)[None, :]
            for months in SEASON_MONTHS
        ]).astype(np.float32)
    seasonal_precip = np.vstack([
        (np.nansum(monthly_precip[list(months)], axis=0) * 4.0)[None, :]
        for months in SEASON_MONTHS
    ]).astype(np.float32)
    valid_precip = np.isfinite(monthly_precip).any(axis=0)
    seasonal_precip[:, ~valid_precip] = np.nan
    annual_precip = np.where(
        valid_precip,
        np.nansum(monthly_precip, axis=0),
        np.nan,
    ).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        annual_temp = np.nanmean(monthly_temp, axis=0).astype(np.float32)
    dry_month_count = np.sum(
        np.where(np.isfinite(monthly_precip), monthly_precip < 30.0, False),
        axis=0,
    ).astype(np.float32)
    dry_month_count[~valid_precip] = np.nan
    return {
        "monthly_temperature_C": monthly_temp,
        "monthly_precip_mm": monthly_precip,
        "seasonal_temperature_C": seasonal_temp,
        "seasonal_precip_mm_yr_equiv": seasonal_precip,
        "annual_temperature_C": annual_temp,
        "annual_precip_mm": annual_precip,
        "dry_month_count": dry_month_count,
    }


def worldclim_metrics(
    grid: SphereGrid,
    fields: dict[str, np.ndarray],
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    annual_temp = np.asarray(fields["annual_temperature_C"], dtype=np.float64)
    annual_precip = np.asarray(fields["annual_precip_mm"], dtype=np.float64)
    seasonal_temp = np.asarray(fields["seasonal_temperature_C"], dtype=np.float64)
    dry_month_count = np.asarray(fields["dry_month_count"], dtype=np.float64)
    land = np.isfinite(annual_temp) & np.isfinite(annual_precip)
    weights = area[land]
    temp_amp = np.nanmax(seasonal_temp[:, land], axis=0) - np.nanmin(
        seasonal_temp[:, land],
        axis=0,
    ) if land.any() else np.asarray([])
    metrics = {
        "worldclim_land_cell_fraction": float(np.mean(land)),
        "worldclim_land_area_fraction": (
            float(area[land].sum() / max(area.sum(), 1.0e-12)) if land.any() else 0.0
        ),
        "land_annual_temperature_mean_C": (
            float(np.average(annual_temp[land], weights=weights)) if land.any() else np.nan
        ),
        "land_annual_temperature_p05_C": (
            float(np.nanpercentile(annual_temp[land], 5)) if land.any() else np.nan
        ),
        "land_annual_temperature_p95_C": (
            float(np.nanpercentile(annual_temp[land], 95)) if land.any() else np.nan
        ),
        "land_seasonal_temperature_amplitude_p50_C": (
            float(np.nanpercentile(temp_amp, 50)) if temp_amp.size else np.nan
        ),
        "land_seasonal_temperature_amplitude_p90_C": (
            float(np.nanpercentile(temp_amp, 90)) if temp_amp.size else np.nan
        ),
        "land_annual_precip_mean_mm": (
            float(np.average(annual_precip[land], weights=weights)) if land.any() else np.nan
        ),
        "land_annual_precip_p50_mm": (
            float(np.nanpercentile(annual_precip[land], 50)) if land.any() else np.nan
        ),
        "land_annual_precip_p90_mm": (
            float(np.nanpercentile(annual_precip[land], 90)) if land.any() else np.nan
        ),
        "land_dry_month_count_p50": (
            float(np.nanpercentile(dry_month_count[land], 50)) if land.any() else np.nan
        ),
        "land_dry_month_count_p90": (
            float(np.nanpercentile(dry_month_count[land], 90)) if land.any() else np.nan
        ),
    }
    return metrics


def noaa_psl_available() -> bool:
    return (
        NOAA_PSL_UWND_10M.exists()
        and NOAA_PSL_VWND_10M.exists()
        and NOAA_PSL_SLP.exists()
    )


def _decode_bytes(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    if hasattr(value, "tobytes"):
        raw = value.tobytes()
        try:
            return raw.decode("utf-8", errors="ignore")
        except Exception:
            return str(value)
    return str(value)


def _netcdf_h5_time_year_month(time_values: np.ndarray, units: str) -> tuple[np.ndarray, np.ndarray]:
    if "hours since 1800-01-01" not in units:
        raise ValueError(f"unsupported NetCDF time units: {units!r}")
    base = datetime(1800, 1, 1)
    years: list[int] = []
    months: list[int] = []
    for value in np.asarray(time_values, dtype=np.float64):
        dt = base + timedelta(hours=float(value))
        years.append(dt.year)
        months.append(dt.month)
    return np.asarray(years, dtype=np.int16), np.asarray(months, dtype=np.int8)


def _load_h5_monthly_climatology(
    path: Path,
    variable: str,
    *,
    start_year: int,
    end_year: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        lat = np.asarray(f["lat"], dtype=np.float32)
        lon = np.asarray(f["lon"], dtype=np.float32)
        time = np.asarray(f["time"], dtype=np.float64)
        units = _decode_bytes(f["time"].attrs.get("units", ""))
        years, months = _netcdf_h5_time_year_month(time, units)
        data = np.asarray(f[variable], dtype=np.float32)
        missing = f[variable].attrs.get("missing_value")
        if missing is not None:
            missing_value = float(np.asarray(missing).ravel()[0])
            data = np.where(data == missing_value, np.nan, data)
        data = np.where(data < -1.0e30, np.nan, data)
    window = (years >= int(start_year)) & (years <= int(end_year))
    if not window.any():
        raise ValueError(f"no NOAA PSL months in {start_year}-{end_year} for {path}")
    clim = []
    for month in range(1, 13):
        mask = window & (months == month)
        if not mask.any():
            raise ValueError(f"missing month {month} in {path}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            clim.append(np.nanmean(data[mask], axis=0)[None, :, :])
    return np.asarray(np.vstack(clim), dtype=np.float32), lat, lon


def _seasonal_mean(monthly: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.asarray([
            np.nanmean(monthly[list(months)], axis=0)
            for months in SEASON_MONTHS
        ], dtype=np.float32)


def sample_latlon_coordinate_raster_to_grid(
    raster: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    grid: SphereGrid,
) -> np.ndarray:
    raster = np.asarray(raster)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    lat_idx = np.abs(lat[:, None] - grid.lat[None, :]).argmin(axis=0)
    lon_grid = (grid.lon + 360.0) % 360.0
    delta_lon = np.abs(((lon[:, None] - lon_grid[None, :] + 180.0) % 360.0) - 180.0)
    lon_idx = delta_lon.argmin(axis=0)
    return raster[lat_idx, lon_idx]


def sample_noaa_psl_to_grid(
    grid: SphereGrid,
    *,
    start_year: int = 1991,
    end_year: int = 2020,
) -> dict[str, np.ndarray]:
    if not noaa_psl_available():
        raise FileNotFoundError("NOAA PSL NCEP wind/pressure files are missing")
    u_monthly, wind_lat, wind_lon = _load_h5_monthly_climatology(
        NOAA_PSL_UWND_10M,
        "uwnd",
        start_year=start_year,
        end_year=end_year,
    )
    v_monthly, wind_lat_v, wind_lon_v = _load_h5_monthly_climatology(
        NOAA_PSL_VWND_10M,
        "vwnd",
        start_year=start_year,
        end_year=end_year,
    )
    slp_monthly, slp_lat, slp_lon = _load_h5_monthly_climatology(
        NOAA_PSL_SLP,
        "slp",
        start_year=start_year,
        end_year=end_year,
    )
    if not (np.allclose(wind_lat, wind_lat_v) and np.allclose(wind_lon, wind_lon_v)):
        raise ValueError("NOAA PSL u/v wind grids do not match")

    u_grid = np.vstack([
        sample_latlon_coordinate_raster_to_grid(u_monthly[i], wind_lat, wind_lon, grid)[None, :]
        for i in range(12)
    ]).astype(np.float32)
    v_grid = np.vstack([
        sample_latlon_coordinate_raster_to_grid(v_monthly[i], wind_lat, wind_lon, grid)[None, :]
        for i in range(12)
    ]).astype(np.float32)
    slp_grid = np.vstack([
        sample_latlon_coordinate_raster_to_grid(slp_monthly[i], slp_lat, slp_lon, grid)[None, :]
        for i in range(12)
    ]).astype(np.float32)
    monthly_wind = np.stack([u_grid, v_grid], axis=2)
    seasonal_wind = np.stack([
        _seasonal_mean(u_grid),
        _seasonal_mean(v_grid),
    ], axis=2)
    seasonal_slp = _seasonal_mean(slp_grid)
    annual_slp = np.nanmean(slp_grid, axis=0).astype(np.float32)
    return {
        "monthly_wind_u10_v10": monthly_wind.astype(np.float32),
        "seasonal_wind_u10_v10": seasonal_wind.astype(np.float32),
        "monthly_slp_hPa": slp_grid.astype(np.float32),
        "seasonal_slp_hPa": seasonal_slp.astype(np.float32),
        "annual_slp_hPa": annual_slp,
        "seasonal_slp_anomaly_hPa": (seasonal_slp - annual_slp[None, :]).astype(np.float32),
    }


def noaa_psl_metrics(fields: dict[str, np.ndarray]) -> dict[str, float]:
    seasonal_wind = np.asarray(fields["seasonal_wind_u10_v10"], dtype=np.float64)
    speed = np.linalg.norm(seasonal_wind, axis=2)
    slp_anom = np.asarray(fields["seasonal_slp_anomaly_hPa"], dtype=np.float64)
    return {
        "seasonal_wind_speed_p50_m_s": float(np.nanpercentile(speed, 50)),
        "seasonal_wind_speed_p90_m_s": float(np.nanpercentile(speed, 90)),
        "seasonal_wind_speed_p99_m_s": float(np.nanpercentile(speed, 99)),
        "seasonal_slp_anomaly_abs_p90_hPa": float(np.nanpercentile(np.abs(slp_anom), 90)),
        "seasonal_slp_anomaly_abs_p99_hPa": float(np.nanpercentile(np.abs(slp_anom), 99)),
    }


def noaa_oisst_available() -> bool:
    return NOAA_OISST_SST_LTM.exists() and NOAA_OISST_ICEC_LTM.exists()


def _read_h5_grid_variable(path: Path, variable: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        lat = np.asarray(f["lat"], dtype=np.float32)
        lon = np.asarray(f["lon"], dtype=np.float32)
        data = np.asarray(f[variable], dtype=np.float32)
        missing = f[variable].attrs.get("missing_value")
        if missing is not None:
            missing_value = float(np.asarray(missing).ravel()[0])
            data = np.where(data == missing_value, np.nan, data)
        fill = f[variable].attrs.get("_FillValue")
        if fill is not None:
            fill_value = float(np.asarray(fill).ravel()[0])
            data = np.where(data == fill_value, np.nan, data)
        valid_range = f[variable].attrs.get("valid_range")
        if valid_range is not None:
            lo, hi = [float(x) for x in np.asarray(valid_range).ravel()[:2]]
            data = np.where((data < lo) | (data > hi), np.nan, data)
        data = np.where(data < -1.0e30, np.nan, data)
    return data.astype(np.float32), lat, lon


def load_noaa_oisst_ltm(
    *,
    sst_path: Path = NOAA_OISST_SST_LTM,
    icec_path: Path = NOAA_OISST_ICEC_LTM,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not sst_path.exists() or not icec_path.exists():
        raise FileNotFoundError(f"NOAA OISST files missing: {sst_path} / {icec_path}")
    sst, lat, lon = _read_h5_grid_variable(sst_path, "sst")
    icec, lat_i, lon_i = _read_h5_grid_variable(icec_path, "icec")
    if not (np.allclose(lat, lat_i) and np.allclose(lon, lon_i)):
        raise ValueError("NOAA OISST SST and ice grids do not match")
    return sst, icec, lat, lon


def sample_noaa_oisst_to_grid(
    grid: SphereGrid,
    *,
    ocean_mask: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    monthly_sst_raster, monthly_ice_raster, lat, lon = load_noaa_oisst_ltm()
    monthly_sst = np.vstack([
        sample_latlon_coordinate_raster_to_grid(
            monthly_sst_raster[i],
            lat,
            lon,
            grid,
        )[None, :]
        for i in range(12)
    ]).astype(np.float32)
    monthly_ice = np.vstack([
        sample_latlon_coordinate_raster_to_grid(
            monthly_ice_raster[i],
            lat,
            lon,
            grid,
        )[None, :]
        for i in range(12)
    ]).astype(np.float32)
    monthly_ice = np.where(np.isfinite(monthly_ice), np.clip(monthly_ice, 0.0, 100.0), np.nan)
    if ocean_mask is not None:
        ocean = np.asarray(ocean_mask, dtype=bool)
        monthly_sst[:, ~ocean] = np.nan
        monthly_ice[:, ~ocean] = np.nan
    seasonal_sst = _seasonal_mean(monthly_sst)
    seasonal_ice = _seasonal_mean(monthly_ice)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        annual_sst = np.nanmean(monthly_sst, axis=0).astype(np.float32)
        annual_ice = np.nanmean(monthly_ice, axis=0).astype(np.float32)
    return {
        "monthly_sst_C": monthly_sst,
        "seasonal_sst_C": seasonal_sst,
        "annual_sst_C": annual_sst,
        "monthly_sea_ice_concentration_pct": monthly_ice.astype(np.float32),
        "seasonal_sea_ice_concentration_pct": seasonal_ice.astype(np.float32),
        "annual_sea_ice_concentration_pct": annual_ice,
    }


def noaa_oisst_metrics(
    grid: SphereGrid,
    fields: dict[str, np.ndarray],
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    annual = np.asarray(fields["annual_sst_C"], dtype=np.float64)
    seasonal = np.asarray(fields["seasonal_sst_C"], dtype=np.float64)
    ice = np.asarray(fields["annual_sea_ice_concentration_pct"], dtype=np.float64)
    valid = np.isfinite(annual)
    tropical = valid & (np.abs(grid.lat) <= 23.5)
    polar = valid & (np.abs(grid.lat) >= 60.0)
    amp = np.nanmax(seasonal[:, valid], axis=0) - np.nanmin(seasonal[:, valid], axis=0) if valid.any() else np.asarray([])
    return {
        "sst_valid_area_fraction": (
            float(area[valid].sum() / max(float(area.sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
        "annual_sst_mean_C": (
            float(np.average(annual[valid], weights=area[valid])) if valid.any() else np.nan
        ),
        "tropical_sst_mean_C": (
            float(np.average(annual[tropical], weights=area[tropical])) if tropical.any() else np.nan
        ),
        "polar_ocean_sst_mean_C": (
            float(np.average(annual[polar], weights=area[polar])) if polar.any() else np.nan
        ),
        "seasonal_sst_amplitude_p90_C": (
            float(np.nanpercentile(amp, 90)) if amp.size else np.nan
        ),
        "annual_sea_ice_area_fraction_gt15pct": (
            float(area[valid & (ice >= 15.0)].sum() / max(float(area[valid].sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
    }


def noaa_aoml_drifter_current_available() -> bool:
    try:
        import lerc  # noqa: F401
    except ImportError:
        return False
    return True


def _arcgis_lerc_tile_path(
    *,
    level: int,
    row: int,
    col: int,
    root: Path = NOAA_AOML_DRIFTER_DIR,
) -> Path:
    return root / f"level{level}" / f"r{row:02d}_c{col:02d}.lerc"


def _download_arcgis_lerc_tile(url: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return path
    request = urllib.request.Request(url, headers={"User-Agent": "aevum-earth-reference/1.0"})
    last_error: BaseException | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                data = response.read()
            break
        except urllib.error.HTTPError:
            raise
        except (
            urllib.error.URLError,
            http.client.RemoteDisconnected,
            http.client.IncompleteRead,
            TimeoutError,
            ConnectionResetError,
        ) as exc:
            last_error = exc
            time.sleep(0.75 * (attempt + 1))
    else:
        raise RuntimeError(f"failed to download ArcGIS tile after retries: {url}: {last_error}")
    if len(data) < 64:
        raise RuntimeError(f"short ArcGIS tile response from {url}")
    path.write_bytes(data)
    return path


def load_noaa_aoml_drifter_current_mosaic(
    *,
    level: int = NOAA_AOML_DRIFTER_LEVEL,
    cache_root: Path = NOAA_AOML_DRIFTER_DIR,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load/cache the NOAA/AOML annual drifter current LERC tile mosaic.

    Returns:
        `(bands, lat, lon)` where bands are speed and direction rasters shaped
        `(2, n_lat, n_lon)` and lat/lon are pixel-center coordinates.
    """
    try:
        import lerc
    except ImportError as exc:
        raise RuntimeError("lerc is required to read ArcGIS current tiles") from exc
    if level not in NOAA_AOML_DRIFTER_RESOLUTIONS:
        raise ValueError(f"unsupported drifter tile level: {level}")
    res = float(NOAA_AOML_DRIFTER_RESOLUTIONS[level])
    xmin, ymin, xmax, ymax = NOAA_AOML_DRIFTER_EXTENT
    origin_x, origin_y = NOAA_AOML_DRIFTER_ORIGIN
    tile_world = ARCGIS_TILE_SIZE * res
    eps = 1.0e-6 * tile_world
    min_col = int(np.floor((xmin - origin_x + eps) / tile_world))
    max_col = int(np.ceil((xmax - origin_x) / tile_world)) - 1
    min_row = int(np.floor((origin_y - ymax + eps) / tile_world))
    max_row = int(np.ceil((origin_y - ymin) / tile_world)) - 1
    min_col = max(0, min_col)
    min_row = max(0, min_row)
    n_cols = (max_col - min_col + 1) * ARCGIS_TILE_SIZE
    n_rows = (max_row - min_row + 1) * ARCGIS_TILE_SIZE
    bands = np.full((2, n_rows, n_cols), np.nan, dtype=np.float32)
    valid = np.zeros((n_rows, n_cols), dtype=bool)
    for row in range(min_row, max_row + 1):
        for col in range(min_col, max_col + 1):
            url = f"{NOAA_AOML_DRIFTER_BASE_URL}/tile/{level}/{row}/{col}"
            path = _arcgis_lerc_tile_path(level=level, row=row, col=col, root=cache_root)
            _download_arcgis_lerc_tile(url, path)
            status, arr, mask = lerc.decode(path.read_bytes())
            if status != 0:
                raise RuntimeError(f"LERC decode failed for {path} with status {status}")
            arr = np.asarray(arr, dtype=np.float32)
            if arr.shape != (2, ARCGIS_TILE_SIZE, ARCGIS_TILE_SIZE):
                raise ValueError(f"unexpected LERC tile shape {arr.shape} for {path}")
            mask = np.asarray(mask, dtype=bool)
            rr = (row - min_row) * ARCGIS_TILE_SIZE
            cc = (col - min_col) * ARCGIS_TILE_SIZE
            tile = np.where(mask[None, :, :], arr, np.nan)
            bands[:, rr:rr + ARCGIS_TILE_SIZE, cc:cc + ARCGIS_TILE_SIZE] = tile
            valid[rr:rr + ARCGIS_TILE_SIZE, cc:cc + ARCGIS_TILE_SIZE] |= mask

    lon = origin_x + (np.arange(n_cols, dtype=np.float64) + min_col * ARCGIS_TILE_SIZE + 0.5) * res
    lat = origin_y - (np.arange(n_rows, dtype=np.float64) + min_row * ARCGIS_TILE_SIZE + 0.5) * res
    col_keep = (lon >= xmin) & (lon <= xmax)
    row_keep = (lat >= ymin) & (lat <= ymax)
    bands = bands[:, row_keep, :][:, :, col_keep]
    valid = valid[row_keep, :][:, col_keep]
    bands[:, ~valid] = np.nan
    lon = lon[col_keep]
    lat = lat[row_keep]
    return bands, lat.astype(np.float32), lon.astype(np.float32)


def sample_noaa_aoml_drifter_current_to_grid(grid: SphereGrid) -> dict[str, np.ndarray]:
    bands, lat, lon = load_noaa_aoml_drifter_current_mosaic()
    speed_raster = np.where(bands[0] >= 0.0, bands[0], np.nan).astype(np.float32)
    direction_raster = bands[1].astype(np.float32)
    theta = np.deg2rad(direction_raster.astype(np.float64))
    # ArcGIS Vector-MagDir uses geographic compass bearings; convert to east/north.
    u_raster = (speed_raster * np.sin(theta)).astype(np.float32)
    v_raster = (speed_raster * np.cos(theta)).astype(np.float32)
    speed = sample_latlon_coordinate_raster_to_grid(speed_raster, lat, lon, grid).astype(np.float32)
    direction = sample_latlon_coordinate_raster_to_grid(direction_raster, lat, lon, grid).astype(np.float32)
    u = sample_latlon_coordinate_raster_to_grid(u_raster, lat, lon, grid).astype(np.float32)
    v = sample_latlon_coordinate_raster_to_grid(v_raster, lat, lon, grid).astype(np.float32)
    current = np.stack([u, v], axis=1).astype(np.float32)
    return {
        "surface_current_u_v": current,
        "annual_surface_current_speed_m_s": speed,
        "annual_surface_current_direction_deg": direction,
    }


def noaa_aoml_drifter_current_metrics(
    grid: SphereGrid,
    fields: dict[str, np.ndarray],
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    speed = np.asarray(fields["annual_surface_current_speed_m_s"], dtype=np.float64)
    valid = np.isfinite(speed) & (speed >= 0.0)
    return {
        "current_valid_area_fraction": (
            float(area[valid].sum() / max(float(area.sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
        "current_speed_p50_m_s": (
            float(np.nanpercentile(speed[valid], 50)) if valid.any() else np.nan
        ),
        "current_speed_p90_m_s": (
            float(np.nanpercentile(speed[valid], 90)) if valid.any() else np.nan
        ),
        "current_speed_p99_m_s": (
            float(np.nanpercentile(speed[valid], 99)) if valid.any() else np.nan
        ),
        "swift_current_area_fraction_gt0_3m_s": (
            float(area[valid & (speed > 0.3)].sum() / max(float(area[valid].sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
    }


def _oscar_climatology_cache_path(
    *,
    level: int = OSCAR_ARCGIS_LEVEL,
    cache_root: Path = OSCAR_ARCGIS_DIR,
) -> Path:
    return cache_root / f"oscar_2001_2020_monthly_climatology_level{level}.npz"


def oscar_monthly_current_available(
    *,
    level: int = OSCAR_ARCGIS_LEVEL,
    cache_root: Path = OSCAR_ARCGIS_DIR,
) -> bool:
    try:
        import lerc  # noqa: F401
    except ImportError:
        return False
    path = _oscar_climatology_cache_path(level=level, cache_root=cache_root)
    return path.exists() and path.stat().st_size > 0


def _download_json(url: str, path: Path, *, refresh: bool = False) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not refresh:
        return json.loads(path.read_text())
    request = urllib.request.Request(url, headers={"User-Agent": "aevum-earth-reference/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read().decode("utf-8")
    data = json.loads(payload)
    path.write_text(json.dumps(data, indent=2, default=_json_default))
    return data


def load_oscar_arcgis_slice_times(
    *,
    cache_root: Path = OSCAR_ARCGIS_DIR,
    refresh: bool = False,
) -> list[int]:
    info = _download_json(
        f"{OSCAR_ARCGIS_SERVICE_URL}/multidimensionalInfo?f=json",
        cache_root / "multidimensional_info.json",
        refresh=refresh,
    )
    variables = info.get("multidimensionalInfo", {}).get("variables", [])
    for variable in variables:
        if variable.get("name") != "Ocean Current":
            continue
        for dimension in variable.get("dimensions", []):
            if dimension.get("name") == "StdTime":
                values = [int(round(float(x))) for x in dimension.get("values", [])]
                if not values:
                    raise ValueError("OSCAR ArcGIS StdTime dimension has no values")
                return values
    raise ValueError("OSCAR ArcGIS multidimensionalInfo missing Ocean Current StdTime")


def _oscar_tile_range(level: int) -> tuple[int, int, int, int, int, int]:
    if level not in OSCAR_ARCGIS_RESOLUTIONS:
        raise ValueError(f"unsupported OSCAR ArcGIS tile level: {level}")
    res = float(OSCAR_ARCGIS_RESOLUTIONS[level])
    xmin, ymin, xmax, ymax = OSCAR_ARCGIS_EXTENT
    origin_x, origin_y = OSCAR_ARCGIS_ORIGIN
    tile_world = ARCGIS_TILE_SIZE * res
    eps = 1.0e-6 * tile_world
    min_col = int(np.floor((xmin - origin_x + eps) / tile_world))
    max_col = int(np.ceil((xmax - origin_x) / tile_world)) - 1
    min_row = int(np.floor((origin_y - ymax + eps) / tile_world))
    max_row = int(np.ceil((origin_y - ymin) / tile_world)) - 1
    min_col = max(0, min_col)
    min_row = max(0, min_row)
    return min_row, max_row, min_col, max_col, max_row - min_row + 1, max_col - min_col + 1


def _oscar_arcgis_tile_path(
    *,
    level: int,
    slice_id: int,
    row: int,
    col: int,
    cache_root: Path = OSCAR_ARCGIS_DIR,
) -> Path:
    return (
        cache_root
        / "tiles"
        / f"level{level}"
        / f"slice{slice_id:03d}"
        / f"r{row:02d}_c{col:02d}.lerc"
    )


def _download_oscar_arcgis_tile(
    *,
    level: int,
    slice_id: int,
    row: int,
    col: int,
    cache_root: Path = OSCAR_ARCGIS_DIR,
) -> Path:
    query = urllib.parse.urlencode({"sliceId": int(slice_id), "blankTile": "false"})
    url = f"{OSCAR_ARCGIS_SERVICE_URL}/tile/{level}/{row}/{col}?{query}"
    path = _oscar_arcgis_tile_path(
        level=level, slice_id=slice_id, row=row, col=col, cache_root=cache_root)
    return _download_arcgis_lerc_tile(url, path)


def _load_oscar_tile(
    *,
    level: int,
    slice_id: int,
    row: int,
    col: int,
    cache_root: Path,
) -> tuple[int, int, np.ndarray, np.ndarray]:
    try:
        import lerc
    except ImportError as exc:
        raise RuntimeError("lerc is required to read OSCAR ArcGIS tiles") from exc
    path = _download_oscar_arcgis_tile(
        level=level,
        slice_id=slice_id,
        row=row,
        col=col,
        cache_root=cache_root,
    )
    status, arr, mask = lerc.decode(path.read_bytes())
    if status != 0:
        raise RuntimeError(f"LERC decode failed for {path} with status {status}")
    arr = np.asarray(arr, dtype=np.float32)
    if arr.shape != (2, ARCGIS_TILE_SIZE, ARCGIS_TILE_SIZE):
        raise ValueError(f"unexpected OSCAR tile shape {arr.shape} for {path}")
    valid = np.asarray(mask, dtype=bool)
    if valid.shape != (ARCGIS_TILE_SIZE, ARCGIS_TILE_SIZE):
        raise ValueError(f"unexpected OSCAR tile mask shape {valid.shape} for {path}")
    valid &= np.isfinite(arr[0]) & np.isfinite(arr[1])
    valid &= (np.abs(arr[0]) < 10.0) & (np.abs(arr[1]) < 10.0)
    return row, col, arr, valid


def _month_index_from_esri_ms(value_ms: int) -> int:
    return datetime.utcfromtimestamp(float(value_ms) / 1000.0).month - 1


def load_oscar_monthly_current_climatology(
    *,
    level: int = OSCAR_ARCGIS_LEVEL,
    cache_root: Path = OSCAR_ARCGIS_DIR,
    allow_download: bool = False,
    max_workers: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Load or build a 12-month OSCAR surface-current climatology.

    The source ArcGIS item exposes 240 monthly U/V slices from 2001-2020.
    This function averages same-calendar-month slices into a compact local
    cache, returning `(u_monthly, v_monthly, y, x, metadata)` where monthly
    arrays have shape `(12, n_y, n_x)` in Web Mercator tile coordinates.
    """
    cache_path = _oscar_climatology_cache_path(level=level, cache_root=cache_root)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with np.load(cache_path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            return (
                data["monthly_u_m_s"].astype(np.float32),
                data["monthly_v_m_s"].astype(np.float32),
                data["y_m"].astype(np.float64),
                data["x_m"].astype(np.float64),
                metadata,
            )
    if not allow_download:
        raise FileNotFoundError(
            f"OSCAR monthly climatology cache missing: {cache_path}. "
            "Run earth-climate-reference with --download-oscar to build it."
        )

    slice_times = load_oscar_arcgis_slice_times(cache_root=cache_root)
    min_row, max_row, min_col, max_col, n_tile_rows, n_tile_cols = _oscar_tile_range(level)
    n_rows = n_tile_rows * ARCGIS_TILE_SIZE
    n_cols = n_tile_cols * ARCGIS_TILE_SIZE
    accum_u = np.zeros((12, n_rows, n_cols), dtype=np.float64)
    accum_v = np.zeros((12, n_rows, n_cols), dtype=np.float64)
    counts = np.zeros((12, n_rows, n_cols), dtype=np.uint16)

    workers = max(1, min(int(max_workers), 5))
    tile_rows = list(range(min_row, max_row + 1))
    tile_cols = list(range(min_col, max_col + 1))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for slice_id, time_ms in enumerate(slice_times):
            month = _month_index_from_esri_ms(time_ms)
            futures = [
                pool.submit(
                    _load_oscar_tile,
                    level=level,
                    slice_id=slice_id,
                    row=row,
                    col=col,
                    cache_root=cache_root,
                )
                for row in tile_rows
                for col in tile_cols
            ]
            for future in futures:
                row, col, arr, valid = future.result()
                rr = (row - min_row) * ARCGIS_TILE_SIZE
                cc = (col - min_col) * ARCGIS_TILE_SIZE
                target = valid
                if not target.any():
                    continue
                accum_u[month, rr:rr + ARCGIS_TILE_SIZE, cc:cc + ARCGIS_TILE_SIZE][target] += arr[0][target]
                accum_v[month, rr:rr + ARCGIS_TILE_SIZE, cc:cc + ARCGIS_TILE_SIZE][target] += arr[1][target]
                counts[month, rr:rr + ARCGIS_TILE_SIZE, cc:cc + ARCGIS_TILE_SIZE][target] += 1

    with np.errstate(invalid="ignore", divide="ignore"):
        monthly_u = (accum_u / counts).astype(np.float32)
        monthly_v = (accum_v / counts).astype(np.float32)
    monthly_u[counts == 0] = np.nan
    monthly_v[counts == 0] = np.nan

    res = float(OSCAR_ARCGIS_RESOLUTIONS[level])
    origin_x, origin_y = OSCAR_ARCGIS_ORIGIN
    x = origin_x + (np.arange(n_cols, dtype=np.float64) + min_col * ARCGIS_TILE_SIZE + 0.5) * res
    y = origin_y - (np.arange(n_rows, dtype=np.float64) + min_row * ARCGIS_TILE_SIZE + 0.5) * res
    metadata = {
        "source_id": "NASA_JPL_OSCAR",
        "source_item": OSCAR_ARCGIS_ITEM_URL,
        "service_url": OSCAR_ARCGIS_SERVICE_URL,
        "level": int(level),
        "slice_count": int(len(slice_times)),
        "tile_rows": [int(min_row), int(max_row)],
        "tile_cols": [int(min_col), int(max_col)],
        "month_sample_counts": [int(np.max(counts[i])) for i in range(12)],
        "cache_path": str(cache_path.relative_to(PROJECT_ROOT)),
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        monthly_u_m_s=monthly_u,
        monthly_v_m_s=monthly_v,
        count=counts,
        y_m=y.astype(np.float32),
        x_m=x.astype(np.float32),
        slice_times_ms=np.asarray(slice_times, dtype=np.int64),
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    return monthly_u, monthly_v, y, x, metadata


def _lonlat_to_webmercator(lon_deg: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = ((np.asarray(lon_deg, dtype=np.float64) + 180.0) % 360.0) - 180.0
    lat = np.clip(
        np.asarray(lat_deg, dtype=np.float64),
        -OSCAR_WEBMERCATOR_MAX_LAT_DEG,
        OSCAR_WEBMERCATOR_MAX_LAT_DEG,
    )
    x = OSCAR_WEBMERCATOR_RADIUS_M * np.deg2rad(lon)
    y = OSCAR_WEBMERCATOR_RADIUS_M * np.log(np.tan(np.pi / 4.0 + np.deg2rad(lat) / 2.0))
    return x, y


def sample_webmercator_raster_to_grid(
    raster: np.ndarray,
    y_m: np.ndarray,
    x_m: np.ndarray,
    grid: SphereGrid,
) -> np.ndarray:
    raster = np.asarray(raster)
    if raster.ndim != 2:
        raise ValueError("raster must be a 2-D Web Mercator array")
    y = np.asarray(y_m, dtype=np.float64)
    x = np.asarray(x_m, dtype=np.float64)
    if raster.shape != (y.size, x.size):
        raise ValueError(f"raster shape {raster.shape} does not match y/x {(y.size, x.size)}")
    dx = float(np.median(np.diff(x)))
    dy = float(abs(np.median(np.diff(y))))
    grid_x, grid_y = _lonlat_to_webmercator(grid.lon, grid.lat)
    col = np.rint((grid_x - float(x[0])) / dx).astype(np.int64)
    row = np.rint((float(y[0]) - grid_y) / dy).astype(np.int64)
    out = np.full(grid.n, np.nan, dtype=np.float32)
    valid = (row >= 0) & (row < y.size) & (col >= 0) & (col < x.size)
    out[valid] = raster[row[valid], col[valid]].astype(np.float32)
    return out


def sample_oscar_monthly_current_to_grid(
    grid: SphereGrid,
    *,
    allow_download: bool = False,
    level: int = OSCAR_ARCGIS_LEVEL,
) -> dict[str, np.ndarray | dict[str, Any]]:
    u_monthly, v_monthly, y, x, metadata = load_oscar_monthly_current_climatology(
        level=level,
        allow_download=allow_download,
    )
    u_cells = np.vstack([
        sample_webmercator_raster_to_grid(u_monthly[i], y, x, grid)[None, :]
        for i in range(12)
    ]).astype(np.float32)
    v_cells = np.vstack([
        sample_webmercator_raster_to_grid(v_monthly[i], y, x, grid)[None, :]
        for i in range(12)
    ]).astype(np.float32)
    monthly_uv = np.stack([u_cells, v_cells], axis=2).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        seasonal_uv = np.stack([
            np.nanmean(monthly_uv[np.asarray(months, dtype=np.int64)], axis=0)
            for months in SEASON_MONTHS
        ], axis=0).astype(np.float32)
        annual_uv = np.nanmean(monthly_uv, axis=0).astype(np.float32)
    monthly_speed = np.linalg.norm(monthly_uv, axis=2).astype(np.float32)
    seasonal_speed = np.linalg.norm(seasonal_uv, axis=2).astype(np.float32)
    annual_speed = np.linalg.norm(annual_uv, axis=1).astype(np.float32)
    return {
        "monthly_surface_current_u_v": monthly_uv,
        "seasonal_surface_current_u_v": seasonal_uv,
        "annual_surface_current_u_v": annual_uv,
        "monthly_surface_current_speed_m_s": monthly_speed,
        "seasonal_surface_current_speed_m_s": seasonal_speed,
        "annual_surface_current_speed_m_s": annual_speed,
        "metadata": metadata,
    }


def oscar_monthly_current_metrics(
    grid: SphereGrid,
    fields: dict[str, Any],
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    annual = np.asarray(fields["annual_surface_current_speed_m_s"], dtype=np.float64)
    monthly = np.asarray(fields["monthly_surface_current_speed_m_s"], dtype=np.float64)
    seasonal = np.asarray(fields["seasonal_surface_current_speed_m_s"], dtype=np.float64)
    valid = np.isfinite(annual) & (annual >= 0.0)
    monthly_valid = np.isfinite(monthly)
    seasonal_amp = (
        np.nanmax(seasonal[:, valid], axis=0) - np.nanmin(seasonal[:, valid], axis=0)
        if valid.any() else np.asarray([])
    )
    return {
        "current_valid_area_fraction": (
            float(area[valid].sum() / max(float(area.sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
        "annual_current_speed_p50_m_s": (
            float(np.nanpercentile(annual[valid], 50)) if valid.any() else np.nan
        ),
        "annual_current_speed_p90_m_s": (
            float(np.nanpercentile(annual[valid], 90)) if valid.any() else np.nan
        ),
        "annual_current_speed_p99_m_s": (
            float(np.nanpercentile(annual[valid], 99)) if valid.any() else np.nan
        ),
        "monthly_current_speed_p90_m_s": (
            float(np.nanpercentile(monthly[monthly_valid], 90)) if monthly_valid.any() else np.nan
        ),
        "seasonal_current_speed_amplitude_p90_m_s": (
            float(np.nanpercentile(seasonal_amp, 90)) if seasonal_amp.size else np.nan
        ),
        "swift_current_area_fraction_gt0_3m_s": (
            float(area[valid & (annual > 0.3)].sum() / max(float(area[valid].sum()), 1.0e-12))
            if valid.any() else 0.0
        ),
    }


def _esa_cci_land_cover_cache_path(
    *,
    cache_root: Path = ESA_CCI_LAND_COVER_DIR,
) -> Path:
    return cache_root / "esa_cci_lc_2020_preview_mosaic.npz"


def esa_cci_land_cover_available(
    *,
    cache_root: Path = ESA_CCI_LAND_COVER_DIR,
) -> bool:
    path = _esa_cci_land_cover_cache_path(cache_root=cache_root)
    return path.exists() and path.stat().st_size > 0


def _hex_to_rgb(value: str) -> tuple[int, int, int]:
    text = value.strip().lstrip("#")
    if len(text) != 6:
        raise ValueError(f"invalid RGB hex color: {value!r}")
    return int(text[0:2], 16), int(text[2:4], 16), int(text[4:6], 16)


def _download_url(url: str, path: Path, *, timeout: int = 120) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return path
    request = urllib.request.Request(url, headers={"User-Agent": "aevum-earth-reference/1.0"})
    last_error: BaseException | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                data = response.read()
            break
        except urllib.error.HTTPError:
            raise
        except (
            urllib.error.URLError,
            http.client.RemoteDisconnected,
            TimeoutError,
            ConnectionResetError,
        ) as exc:
            last_error = exc
            time.sleep(0.75 * (attempt + 1))
    else:
        raise RuntimeError(f"failed to download URL after retries: {url}: {last_error}")
    if len(data) < 64:
        raise RuntimeError(f"short response from {url}")
    path.write_bytes(data)
    return path


def load_esa_cci_land_cover_items(
    *,
    year: int = ESA_CCI_LAND_COVER_YEAR,
    cache_root: Path = ESA_CCI_LAND_COVER_DIR,
    refresh: bool = False,
) -> list[dict[str, Any]]:
    payload = _download_json(
        ESA_CCI_LAND_COVER_ITEMS_URL,
        cache_root / f"esa_cci_lc_items_{year}.json",
        refresh=refresh,
    )
    items = [
        item for item in payload.get("features", [])
        if f"-P1Y-{int(year)}-" in str(item.get("id", ""))
    ]
    if len(items) != 32:
        raise ValueError(f"expected 32 ESA CCI LC {year} preview tiles, found {len(items)}")
    return items


def _esa_cci_lccs_to_broad_class(value: int) -> int:
    code = int(value)
    if code == 210:
        return 1
    if code in {10, 11, 12, 20, 30}:
        return 2
    if code in {50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 160, 170}:
        return 3
    if code in {40, 100, 110, 120, 121, 122, 130, 140}:
        return 4
    if code == 180:
        return 5
    if code == 190:
        return 6
    if code in {150, 151, 152, 153, 200, 201, 202}:
        return 7
    if code == 220:
        return 8
    return 0


def _esa_cci_color_to_lccs_map(items: list[dict[str, Any]]) -> dict[tuple[int, int, int], int]:
    classes = items[0]["assets"]["lccs_class"].get("classification:classes", [])
    color_map: dict[tuple[int, int, int], int] = {}
    for entry in classes:
        color = entry.get("color_hint") or entry.get("color-hint")
        if not color:
            continue
        rgb = _hex_to_rgb(str(color))
        value = int(entry.get("value", 0))
        if rgb not in color_map or value < color_map[rgb]:
            color_map[rgb] = value
    return color_map


def load_esa_cci_land_cover_preview_mosaic(
    *,
    year: int = ESA_CCI_LAND_COVER_YEAR,
    cache_root: Path = ESA_CCI_LAND_COVER_DIR,
    allow_download: bool = False,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Load or build a coarse global ESA CCI/C3S land-cover mosaic.

    Planetary Computer hosts 45-degree COG tiles and provides rendered
    categorical previews.  This adapter uses those previews plus STAC class
    color hints to produce a compact, same-grid calibration reference.  It is
    a coarse land-cover envelope, not a substitute for the original 300 m data.
    """
    cache_path = _esa_cci_land_cover_cache_path(cache_root=cache_root)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with np.load(cache_path, allow_pickle=False) as data:
            return (
                data["lccs_class"].astype(np.uint16),
                data["broad_class"].astype(np.uint8),
                json.loads(str(data["metadata_json"])),
            )
    if not allow_download:
        raise FileNotFoundError(
            f"ESA CCI land-cover preview cache missing: {cache_path}. "
            "Run earth-climate-reference with --download-land-cover to build it."
        )

    items = load_esa_cci_land_cover_items(year=year, cache_root=cache_root)
    color_to_lccs = _esa_cci_color_to_lccs_map(items)
    color_lookup = np.zeros(1 << 24, dtype=np.uint16)
    for rgb, value in color_to_lccs.items():
        key = (int(rgb[0]) << 16) | (int(rgb[1]) << 8) | int(rgb[2])
        color_lookup[key] = int(value)

    tile = ESA_CCI_LAND_COVER_PREVIEW_SIZE
    rows = int(round(180.0 / ESA_CCI_LAND_COVER_TILE_DEG))
    cols = int(round(360.0 / ESA_CCI_LAND_COVER_TILE_DEG))
    lccs = np.zeros((rows * tile, cols * tile), dtype=np.uint16)
    for item in items:
        item_id = str(item["id"])
        href = item["assets"]["rendered_preview"]["href"]
        path = cache_root / "previews" / f"{item_id}.png"
        _download_url(href, path)
        with Image.open(path) as image:
            rgb = np.asarray(image.convert("RGBA"), dtype=np.uint8)
        if rgb.shape[0] != tile or rgb.shape[1] != tile:
            raise ValueError(f"unexpected ESA CCI preview size {rgb.shape} for {item_id}")
        key = (
            (rgb[:, :, 0].astype(np.uint32) << 16)
            | (rgb[:, :, 1].astype(np.uint32) << 8)
            | rgb[:, :, 2].astype(np.uint32)
        )
        arr = color_lookup[key].astype(np.uint16)
        arr[rgb[:, :, 3] == 0] = 0
        lon_min, lat_min, lon_max, lat_max = [float(x) for x in item["bbox"]]
        rr = int(round((90.0 - lat_max) / ESA_CCI_LAND_COVER_TILE_DEG)) * tile
        cc = int(round((lon_min + 180.0) / ESA_CCI_LAND_COVER_TILE_DEG)) * tile
        lccs[rr:rr + tile, cc:cc + tile] = arr

    broad_lookup = np.zeros(256, dtype=np.uint8)
    for value in range(broad_lookup.size):
        broad_lookup[value] = _esa_cci_lccs_to_broad_class(value)
    broad = broad_lookup[np.clip(lccs, 0, 255)].astype(np.uint8)
    metadata = {
        "source_id": "ESA_CCI_LAND_COVER",
        "source_collection": ESA_CCI_LAND_COVER_COLLECTION_URL,
        "year": int(year),
        "tile_count": int(len(items)),
        "preview_size": int(tile),
        "mosaic_shape": [int(lccs.shape[0]), int(lccs.shape[1])],
        "broad_class_names": LAND_COVER_BROAD_NAMES,
        "cache_path": str(cache_path.relative_to(PROJECT_ROOT)),
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        lccs_class=lccs,
        broad_class=broad,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True, default=_json_default)),
    )
    return lccs, broad, metadata


def sample_esa_cci_land_cover_to_grid(
    grid: SphereGrid,
    *,
    allow_download: bool = False,
) -> dict[str, Any]:
    lccs, broad, metadata = load_esa_cci_land_cover_preview_mosaic(
        allow_download=allow_download)
    lccs_cells = sample_latlon_raster_to_grid(
        lccs, grid, north_to_south=True, lon_origin_deg=-180.0).astype(np.uint16)
    broad_cells = sample_latlon_raster_to_grid(
        broad, grid, north_to_south=True, lon_origin_deg=-180.0).astype(np.uint8)
    return {
        "esa_cci_lccs_class": lccs_cells,
        "esa_cci_land_cover_broad_class": broad_cells,
        "metadata": metadata,
    }


def esa_cci_land_cover_metrics(
    grid: SphereGrid,
    fields: dict[str, Any],
    *,
    elevation_m: np.ndarray | None = None,
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    lccs = np.asarray(fields["esa_cci_lccs_class"], dtype=np.uint16)
    broad = np.asarray(fields["esa_cci_land_cover_broad_class"], dtype=np.uint8)
    valid = lccs > 0
    metrics: dict[str, float] = {
        "valid_area_fraction": float(area[valid].sum() / total_area),
        "observed_lccs_class_count": float(len(set(int(x) for x in lccs[valid]))),
    }
    for code, name in LAND_COVER_BROAD_NAMES.items():
        if code == 0:
            continue
        key = name.replace("-", "_").replace(" ", "_")
        metrics[f"{key}_area_fraction"] = float(area[broad == int(code)].sum() / total_area)
    if elevation_m is not None:
        elevation = np.asarray(elevation_m, dtype=np.float64)
        water = broad == 1
        ocean = elevation < 0.0
        metrics["water_ocean_class_mismatch_area_fraction"] = float(
            area[valid & (water != ocean)].sum() / total_area
        )
    return metrics


def koppen_available() -> bool:
    return gloh2o_koppen_available() or KOPPEN_ASCII_ZIP.exists()


def gloh2o_koppen_available() -> bool:
    return GLOH2O_KOPPEN_ZIP.exists()


def load_gloh2o_koppen_raster(
    *,
    path: Path = GLOH2O_KOPPEN_ZIP,
    period: str = GLOH2O_KOPPEN_PERIOD,
    resolution: str = GLOH2O_KOPPEN_RESOLUTION,
) -> tuple[np.ndarray, dict[int, str]]:
    """Load Beck et al. / GloH2O Koppen-Geiger GeoTIFF from the cached zip."""
    if not path.exists():
        raise FileNotFoundError(path)
    member = f"{period}/koppen_geiger_{resolution}.tif"
    with zipfile.ZipFile(path) as zf:
        if member not in zf.namelist():
            raise ValueError(f"missing {member} in {path}")
        data = zf.read(member)
    with Image.open(BytesIO(data)) as image:
        raster = np.asarray(image).astype(np.int16)
    return raster, dict(GLOH2O_KOPPEN_CODE_TO_CLASS)


def load_koppen_ascii(
    path: Path = KOPPEN_ASCII_ZIP,
) -> tuple[np.ndarray, dict[int, str]]:
    """Load the 0.5-degree Kottek/Rubel Koppen-Geiger ASCII grid."""
    if not path.exists():
        raise FileNotFoundError(path)
    raster = np.zeros((360, 720), dtype=np.int16)
    with zipfile.ZipFile(path) as zf:
        name = next(n for n in zf.namelist() if n.endswith(".txt"))
        with zf.open(name) as f:
            header = f.readline()
            if b"Cls" not in header:
                raise ValueError(f"unexpected Koppen ASCII header: {header!r}")
            for line in f:
                parts = line.decode("ascii", errors="ignore").split()
                if len(parts) != 3:
                    continue
                lat = float(parts[0])
                lon = float(parts[1])
                cls = parts[2]
                row = int(round((lat + 89.75) / 0.5))
                col = int(round((lon + 179.75) / 0.5))
                if 0 <= row < 360 and 0 <= col < 720:
                    raster[row, col] = KOPPEN_TO_CODE.get(cls, 0)
    return raster, dict(KOPPEN_CODE_TO_CLASS)


def _koppen_major_from_code(
    code: np.ndarray,
    class_names: dict[int, str] | None = None,
) -> np.ndarray:
    code = np.asarray(code, dtype=np.int16)
    class_names = class_names or KOPPEN_CODE_TO_CLASS
    major = np.zeros_like(code, dtype=np.int16)
    for numeric, label in class_names.items():
        first = label[0]
        if first == "A":
            major[code == numeric] = 1
        elif first == "B":
            major[code == numeric] = 2
        elif first == "C":
            major[code == numeric] = 3
        elif first == "D":
            major[code == numeric] = 4
        elif first == "E":
            major[code == numeric] = 5
    return major


def _biome_proxy_from_koppen_code(
    code: np.ndarray,
    class_names: dict[int, str] | None = None,
) -> np.ndarray:
    """Map Koppen classes to Aevum's coarse biome codes for calibration only."""
    code = np.asarray(code, dtype=np.int16)
    class_names = class_names or KOPPEN_CODE_TO_CLASS
    biome = np.zeros_like(code, dtype=np.int16)
    for numeric, label in class_names.items():
        mask = code == numeric
        if label.startswith("A"):
            biome[mask] = 6  # tropical
        elif label.startswith("BW"):
            biome[mask] = 2  # desert
        elif label.startswith("BS"):
            biome[mask] = 3  # grassland / steppe
        elif label.startswith(("C", "D")):
            biome[mask] = 4  # forest / temperate-boreal
        elif label == "ET":
            biome[mask] = 5  # tundra
        elif label == "EF":
            biome[mask] = 1  # ice
    return biome


def sample_koppen_to_grid(
    grid: SphereGrid,
    *,
    elevation_m: np.ndarray | None = None,
    prefer_gloh2o: bool = True,
) -> dict[str, np.ndarray | dict[int, str]]:
    if prefer_gloh2o and gloh2o_koppen_available():
        raster, class_names = load_gloh2o_koppen_raster()
        koppen = sample_latlon_raster_to_grid(
            raster,
            grid,
            north_to_south=True,
            lon_origin_deg=-180.0,
        ).astype(np.int16)
        source_id = "GLOH2O_KOPPEN_GEIGER"
        native_resolution = GLOH2O_KOPPEN_RESOLUTION
        period = GLOH2O_KOPPEN_PERIOD
    else:
        raster, class_names = load_koppen_ascii()
        koppen = sample_latlon_raster_to_grid(
            raster,
            grid,
            north_to_south=False,
            lon_origin_deg=-180.0,
        ).astype(np.int16)
        source_id = "KOTTEK_RUBEL_KOPPEN_2006_ASCII"
        native_resolution = "0p5"
        period = "1951_2000"
    if elevation_m is not None:
        koppen = np.where(np.asarray(elevation_m) >= 0.0, koppen, 0).astype(np.int16)
    major = _koppen_major_from_code(koppen, class_names)
    biome_proxy = _biome_proxy_from_koppen_code(koppen, class_names)
    return {
        "koppen_class": koppen,
        "koppen_major_class": major,
        "biome_class_proxy": biome_proxy,
        "koppen_class_names": class_names,
        "koppen_source_id": source_id,
        "koppen_native_resolution": native_resolution,
        "koppen_period": period,
    }


def resolve_ecoregions_available() -> bool:
    return RESOLVE_ECOREGIONS_ZIP.exists()


def sample_resolve_ecoregions_to_grid(
    grid: SphereGrid,
    *,
    land_mask: np.ndarray | None = None,
    path: Path = RESOLVE_ECOREGIONS_ZIP,
) -> dict[str, Any]:
    """Sample RESOLVE Ecoregions 2017 polygons to Aevum cell centers."""
    if not path.exists():
        raise FileNotFoundError(path)
    try:
        import shapefile  # type: ignore
    except ImportError as exc:
        raise RuntimeError("pyshp is required to read RESOLVE shapefiles") from exc

    lon = ((np.asarray(grid.lon, dtype=np.float64) + 180.0) % 360.0) - 180.0
    lat = np.asarray(grid.lat, dtype=np.float64)
    points = np.column_stack([lon, lat])
    active = np.ones(grid.n, dtype=bool)
    if land_mask is not None:
        active &= np.asarray(land_mask, dtype=bool)
    biome = np.zeros(grid.n, dtype=np.int16)
    ecoregion = np.zeros(grid.n, dtype=np.int32)
    ecoregion_names: dict[int, str] = {}
    biome_names: dict[int, str] = dict(RESOLVE_BIOME_NAMES)

    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(path) as zf:
            zf.extractall(tmp)
        reader = shapefile.Reader(str(Path(tmp) / "Ecoregions2017.shp"), encoding="latin1")
        for shape_record in reader.iterShapeRecords():
            record = shape_record.record.as_dict()
            biome_id = int(record["BIOME_NUM"])
            eco_id = int(record["ECO_ID"])
            if eco_id <= 0:
                eco_id = int(record["OBJECTID"])
            ecoregion_names[eco_id] = str(record["ECO_NAME"])
            biome_names[biome_id] = str(record["BIOME_NAME"])
            shape = shape_record.shape
            xmin, ymin, xmax, ymax = [float(x) for x in shape.bbox]
            candidate_base = (
                active
                & (lon >= xmin - 1.0e-9)
                & (lon <= xmax + 1.0e-9)
                & (lat >= ymin - 1.0e-9)
                & (lat <= ymax + 1.0e-9)
            )
            if not candidate_base.any():
                continue
            part_starts = list(shape.parts) + [len(shape.points)]
            assigned = np.zeros(grid.n, dtype=bool)
            for start, stop in zip(part_starts[:-1], part_starts[1:]):
                vertices = np.asarray(shape.points[start:stop], dtype=np.float64)
                if vertices.shape[0] < 3:
                    continue
                pxmin, pymin = np.nanmin(vertices, axis=0)
                pxmax, pymax = np.nanmax(vertices, axis=0)
                candidate = (
                    candidate_base
                    & (lon >= pxmin - 1.0e-9)
                    & (lon <= pxmax + 1.0e-9)
                    & (lat >= pymin - 1.0e-9)
                    & (lat <= pymax + 1.0e-9)
                )
                candidate_idx = np.flatnonzero(candidate)
                if candidate_idx.size == 0:
                    continue
                inside = MplPath(vertices).contains_points(points[candidate_idx])
                if inside.any():
                    assigned[candidate_idx[inside]] = True
            if assigned.any():
                biome[assigned] = biome_id
                ecoregion[assigned] = eco_id
                active[assigned] = False

    return {
        "resolve_biome_class": biome,
        "resolve_ecoregion_id": ecoregion,
        "resolve_biome_names": biome_names,
        "resolve_ecoregion_names": ecoregion_names,
    }


def categorical_fraction_metrics(
    grid: SphereGrid,
    field: np.ndarray,
    *,
    prefix: str,
) -> dict[str, float]:
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total = max(float(area.sum()), 1.0e-12)
    out: dict[str, float] = {}
    for value in sorted(int(x) for x in np.unique(field)):
        out[f"{prefix}_{value}_area_fraction"] = float(area[field == value].sum() / total)
    return out


def resolve_ecoregion_metrics(
    grid: SphereGrid,
    fields: dict[str, Any],
    *,
    land_mask: np.ndarray | None = None,
) -> dict[str, float]:
    biome = np.asarray(fields["resolve_biome_class"], dtype=np.int16)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    if land_mask is None:
        land = np.ones(grid.n, dtype=bool)
    else:
        land = np.asarray(land_mask, dtype=bool)
    land_area = max(float(area[land].sum()), 1.0e-12)
    assigned = biome > 0
    metrics = categorical_fraction_metrics(grid, biome, prefix="resolve_biome")
    metrics["resolve_assigned_land_area_fraction"] = float(
        area[land & assigned].sum() / land_area
    )
    metrics["resolve_assigned_cell_fraction"] = float(np.mean(assigned))
    metrics["resolve_observed_biome_count"] = float(len(set(int(x) for x in biome if x > 0)))
    return metrics



def earth_elevation_metrics(
    grid: SphereGrid,
    elevation_m: np.ndarray,
) -> dict[str, float]:
    elevation_m = np.asarray(elevation_m, dtype=np.float64)
    area = np.asarray(grid.cell_area, dtype=np.float64)
    total_area = max(float(area.sum()), 1.0e-12)
    land = elevation_m >= 0.0
    ocean = ~land
    land_area = max(float(area[land].sum()), 1.0e-12)
    ocean_area = max(float(area[ocean].sum()), 1.0e-12)
    depth = np.where(ocean, -elevation_m, 0.0)
    return {
        "land_fraction": float(area[land].sum() / total_area),
        "ocean_fraction": float(area[ocean].sum() / total_area),
        "land_elevation_mean_m": (
            float(np.average(elevation_m[land], weights=area[land])) if land.any() else 0.0
        ),
        "land_elevation_p50_m": float(np.percentile(elevation_m[land], 50)) if land.any() else 0.0,
        "land_elevation_p95_m": float(np.percentile(elevation_m[land], 95)) if land.any() else 0.0,
        "high_land_fraction_gt2500m": float(area[land & (elevation_m > 2500.0)].sum() / land_area),
        "lowland_fraction_lt500m": float(area[land & (elevation_m < 500.0)].sum() / land_area),
        "lowland_fraction_lt1000m": float(area[land & (elevation_m < 1000.0)].sum() / land_area),
        "ocean_depth_mean_m": (
            float(np.average(depth[ocean], weights=area[ocean])) if ocean.any() else 0.0
        ),
        "ocean_depth_p50_m": float(np.percentile(depth[ocean], 50)) if ocean.any() else 0.0,
        "ocean_depth_p95_m": float(np.percentile(depth[ocean], 95)) if ocean.any() else 0.0,
        "shelf_fraction_of_ocean": float(area[ocean & (depth <= 200.0)].sum() / ocean_area),
        "slope_rise_fraction_of_ocean": float(
            area[ocean & (depth > 200.0) & (depth <= 3500.0)].sum() / ocean_area
        ),
        "abyss_fraction_of_ocean": float(
            area[ocean & (depth > 3500.0) & (depth <= 6000.0)].sum() / ocean_area
        ),
        "trench_and_hadal_fraction_of_ocean": float(
            area[ocean & (depth > 6000.0)].sum() / ocean_area
        ),
    }


def _render_elevation_map(
    grid: SphereGrid,
    elevation_m: np.ndarray,
    out_path: Path,
    *,
    title: str,
    width: int,
    height: int,
) -> Path:
    raster = render.to_raster_continuous(
        grid,
        np.asarray(elevation_m, dtype=np.float64),
        width=int(width),
        height=int(height),
        preserve_sign=True,
    )
    fig, ax = plt.subplots(figsize=(10, 5))
    im = render.render_elevation_raster(ax, raster, title=title)
    render.add_elevation_colorbar(fig, ax, im)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_scalar_map(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap,
    width: int,
    height: int,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    raster = render.to_raster(grid, np.asarray(field), width=width, height=height)
    masked = np.ma.masked_invalid(raster)
    if vmin is None or vmax is None:
        values = np.asarray(field)[np.isfinite(field)]
        if values.size:
            if vmin is None:
                vmin = float(np.nanpercentile(values, 2))
            if vmax is None:
                vmax = float(np.nanpercentile(values, 98))
        else:
            vmin = 0.0 if vmin is None else vmin
            vmax = 1.0 if vmax is None else vmax
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(
        masked,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.72)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_seasonal_panels(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap,
    width: int,
    height: int,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    field = np.asarray(field)
    rasters = [
        render.to_raster(grid, field[i], width=width, height=height)
        for i in range(4)
    ]
    combined = np.concatenate([r[np.isfinite(r)].ravel() for r in rasters])
    if combined.size:
        if vmin is None:
            vmin = float(np.nanpercentile(combined, 2))
        if vmax is None:
            vmax = float(np.nanpercentile(combined, 98))
    else:
        vmin = 0.0 if vmin is None else vmin
        vmax = 1.0 if vmax is None else vmax
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), constrained_layout=True)
    im = None
    for ax, label, raster in zip(axes.ravel(), SEASON_LABELS, rasters):
        im = ax.imshow(
            np.ma.masked_invalid(raster),
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            extent=[-180, 180, -90, 90],
        )
        ax.set_title(f"{label} {title}")
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_koppen_major_map(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    width: int,
    height: int,
) -> Path:
    raster = render.to_raster(grid, np.asarray(field), width=width, height=height)
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(
        raster,
        cmap=KOPPEN_MAJOR_CMAP,
        norm=KOPPEN_MAJOR_NORM,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, shrink=0.72, ticks=[0, 1, 2, 3, 4, 5])
    cb.ax.set_yticklabels(["ocean", "A", "B", "C", "D", "E"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_koppen_class_map(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    width: int,
    height: int,
) -> Path:
    raster = render.to_raster(grid, np.asarray(field), width=width, height=height)
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(
        raster,
        cmap=GLOH2O_KOPPEN_CMAP,
        norm=GLOH2O_KOPPEN_NORM,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(title)
    ticks = [0, 1, 4, 8, 14, 17, 25, 29, 30]
    cb = fig.colorbar(im, ax=ax, shrink=0.72, ticks=ticks)
    cb.ax.set_yticklabels(["ocean", "Af", "BWh", "Csa", "Cfa", "Dsa", "Dfa", "ET", "EF"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_biome_proxy_map(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    width: int,
    height: int,
) -> Path:
    raster = render.to_raster(grid, np.asarray(field), width=width, height=height)
    cmap = ListedColormap(render.BIOME_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(render.BIOME_COLORS) + 0.5, 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(raster, cmap=cmap, norm=norm, extent=[-180, 180, -90, 90])
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, shrink=0.72, ticks=list(range(len(render.BIOME_COLORS))))
    cb.ax.set_yticklabels(["ocean", "ice", "desert", "grass", "forest", "tundra", "tropical"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_resolve_biome_map(
    grid: SphereGrid,
    field: np.ndarray,
    out_path: Path,
    *,
    title: str,
    width: int,
    height: int,
) -> Path:
    raster = render.to_raster(grid, np.asarray(field), width=width, height=height)
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(
        raster,
        cmap=RESOLVE_BIOME_CMAP,
        norm=RESOLVE_BIOME_NORM,
        extent=[-180, 180, -90, 90],
    )
    ax.set_title(title)
    ticks = list(range(len(RESOLVE_BIOME_COLORS)))
    cb = fig.colorbar(im, ax=ax, shrink=0.72, ticks=ticks)
    cb.ax.set_yticklabels([str(tick) for tick in ticks])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _render_contact_sheet(
    panels: list[tuple[str, Path]],
    out_path: Path,
    *,
    columns: int = 2,
    tile_w: int = 440,
    tile_h: int = 250,
) -> Path | None:
    existing = [(label, path) for label, path in panels if path.exists()]
    if not existing:
        return None
    font = ImageFont.load_default()
    pad = 10
    label_h = 24
    rows = int(np.ceil(len(existing) / max(columns, 1)))
    canvas = Image.new(
        "RGB",
        (
            columns * (tile_w + pad) + pad,
            rows * (tile_h + label_h + pad) + pad,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    for idx, (label, path) in enumerate(existing):
        row = idx // columns
        col = idx % columns
        x = pad + col * (tile_w + pad)
        y = pad + row * (tile_h + label_h + pad)
        draw.text((x + 4, y + 4), label, fill=(0, 0, 0), font=font)
        with Image.open(path) as image:
            im = image.convert("RGB")
            im.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
            bg = Image.new("RGB", (tile_w, tile_h), "white")
            bg.paste(im, ((tile_w - im.width) // 2, (tile_h - im.height) // 2))
            canvas.paste(bg, (x, y + label_h))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)
    return out_path


def _write_arrays(
    out_path: Path,
    grid: SphereGrid,
    elevation_m: np.ndarray,
    metrics: dict[str, float],
    extra_arrays: dict[str, np.ndarray] | None = None,
    extra_metadata: dict[str, Any] | None = None,
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {
        "lat": grid.lat,
        "lon": grid.lon,
        "cell_area": grid.cell_area,
        "earth__elevation_m": np.asarray(elevation_m, dtype=np.float32),
        "earth__land_mask": np.asarray(elevation_m >= 0.0, dtype=bool),
        "metrics_json": np.asarray(json.dumps(metrics, sort_keys=True)),
    }
    for key, value in (extra_arrays or {}).items():
        arrays[key.replace(".", "__")] = np.asarray(value)
    if extra_metadata:
        arrays["metadata_json"] = np.asarray(
            json.dumps(extra_metadata, sort_keys=True, default=_json_default)
        )
    np.savez_compressed(
        out_path,
        **arrays,
    )
    return out_path


def run_earth_climate_reference(
    config: EarthClimateReferenceConfig,
    outdir: Path,
) -> dict[str, Any]:
    """Build the current Earth reference manifest and same-grid assets."""
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = build_source_manifest()
    manifest_path = outdir / "source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=_json_default))
    if config.manifest_out is not None:
        write_source_manifest(Path(config.manifest_out))

    entries: list[dict[str, Any]] = []
    warnings: list[str] = []
    for cells in config.cells:
        grid = SphereGrid.fibonacci(int(cells), CONSTANTS.EARTH_RADIUS)
        elevation = sample_etopo5_to_grid(grid)
        metrics: dict[str, Any] = {
            "elevation": earth_elevation_metrics(grid, elevation)
        }
        extra_arrays: dict[str, np.ndarray] = {}
        metadata: dict[str, Any] = {
            "source_ids": ["NOAA_ETOPO5_LOCAL"],
            "koppen_class_names": KOPPEN_CODE_TO_CLASS,
        }
        asset_paths: list[str] = []
        panel_paths: list[tuple[str, Path]] = []
        if config.render_assets:
            png = _render_elevation_map(
                grid,
                elevation,
                outdir / f"earth_reference_{int(cells)}cells_elevation.png",
                title=f"Earth ETOPO5 sampled to {int(cells)} Aevum cells",
                width=int(config.width),
                height=int(config.height),
            )
            asset_paths.append(str(png))
            panel_paths.append(("Elevation", png))

        if config.include_etopo2022:
            if etopo2022_available():
                etopo2022 = sample_etopo2022_to_grid(grid)
                extra_arrays.update({
                    "earth.etopo2022_elevation_m": etopo2022,
                    "earth.etopo2022_land_mask": etopo2022 >= 0.0,
                    "earth.etopo2022_minus_etopo5_m": etopo2022 - elevation,
                })
                metrics["etopo2022_elevation"] = earth_elevation_metrics(
                    grid, etopo2022)
                metrics["etopo2022_crosscheck"] = etopo2022_crosscheck_metrics(
                    grid, elevation, etopo2022)
                metadata["source_ids"].append("NOAA_ETOPO_2022")
                metadata["etopo2022_cache_path"] = str(
                    ETOPO2022_OPENDAP_ASC.relative_to(PROJECT_ROOT))
                metadata["etopo2022_cache_url"] = ETOPO2022_OPENDAP_ASC_URL
                if config.render_assets:
                    etopo2022_png = _render_elevation_map(
                        grid,
                        etopo2022,
                        outdir / f"earth_reference_{int(cells)}cells_etopo2022_elevation.png",
                        title=f"Earth ETOPO 2022 sampled to {int(cells)} Aevum cells",
                        width=int(config.width),
                        height=int(config.height),
                    )
                    delta_png = _render_scalar_map(
                        grid,
                        etopo2022 - elevation,
                        outdir / f"earth_reference_{int(cells)}cells_etopo2022_minus_etopo5.png",
                        title="ETOPO 2022 minus ETOPO5 elevation (m)",
                        cmap="coolwarm",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-1500.0,
                        vmax=1500.0,
                    )
                    for label, path in (
                        ("ETOPO 2022", etopo2022_png),
                        ("ETOPO 2022 - ETOPO5", delta_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))
            else:
                warnings.append(
                    "ETOPO 2022 OPeNDAP ASCII cache missing; ETOPO5 cross-check skipped"
                )

        if config.include_worldclim:
            if worldclim_available():
                wc = sample_worldclim_to_grid(grid)
                extra_arrays.update({
                    "earth.monthly_temperature_C": wc["monthly_temperature_C"],
                    "earth.monthly_precip_mm": wc["monthly_precip_mm"],
                    "earth.seasonal_temperature_C": wc["seasonal_temperature_C"],
                    "earth.seasonal_precip_mm_yr_equiv": wc[
                        "seasonal_precip_mm_yr_equiv"
                    ],
                    "earth.annual_temperature_C": wc["annual_temperature_C"],
                    "earth.annual_precip_mm": wc["annual_precip_mm"],
                    "earth.dry_month_count": wc["dry_month_count"],
                })
                metrics["worldclim"] = worldclim_metrics(grid, wc)
                metadata["source_ids"].append("WORLDCLIM_2_1")
                if config.render_assets:
                    temp_png = _render_scalar_map(
                        grid,
                        wc["annual_temperature_C"],
                        outdir / f"earth_reference_{int(cells)}cells_temperature.png",
                        title="WorldClim annual land temperature (C)",
                        cmap=render.TEMPERATURE_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-30.0,
                        vmax=32.0,
                    )
                    precip_png = _render_scalar_map(
                        grid,
                        wc["annual_precip_mm"],
                        outdir / f"earth_reference_{int(cells)}cells_precip.png",
                        title="WorldClim annual land precipitation (mm/yr)",
                        cmap=render.PRECIP_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=3000.0,
                    )
                    temp_seasons_png = _render_seasonal_panels(
                        grid,
                        wc["seasonal_temperature_C"],
                        outdir / f"earth_reference_{int(cells)}cells_temperature_seasons.png",
                        title="WorldClim temperature (C)",
                        cmap=render.TEMPERATURE_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-35.0,
                        vmax=35.0,
                    )
                    precip_seasons_png = _render_seasonal_panels(
                        grid,
                        wc["seasonal_precip_mm_yr_equiv"],
                        outdir / f"earth_reference_{int(cells)}cells_precip_seasons.png",
                        title="WorldClim precipitation (mm/yr equiv.)",
                        cmap=render.PRECIP_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=4000.0,
                    )
                    for label, path in (
                        ("Annual Temp", temp_png),
                        ("Annual Precip", precip_png),
                        ("Seasonal Temp", temp_seasons_png),
                        ("Seasonal Precip", precip_seasons_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))
            else:
                warnings.append(
                    "WorldClim zips missing; temperature/precipitation reference skipped"
                )

        if config.include_sst:
            if noaa_oisst_available():
                sst = sample_noaa_oisst_to_grid(grid, ocean_mask=elevation < 0.0)
                extra_arrays.update({
                    "earth.monthly_sst_C": sst["monthly_sst_C"],
                    "earth.seasonal_sst_C": sst["seasonal_sst_C"],
                    "earth.annual_sst_C": sst["annual_sst_C"],
                    "earth.monthly_sea_ice_concentration_pct": sst[
                        "monthly_sea_ice_concentration_pct"
                    ],
                    "earth.seasonal_sea_ice_concentration_pct": sst[
                        "seasonal_sea_ice_concentration_pct"
                    ],
                    "earth.annual_sea_ice_concentration_pct": sst[
                        "annual_sea_ice_concentration_pct"
                    ],
                })
                metrics["noaa_oisst_v2"] = noaa_oisst_metrics(grid, sst)
                metadata["source_ids"].append("NOAA_OISST_V2_LTM_1991_2020")
                if config.render_assets:
                    sst_png = _render_scalar_map(
                        grid,
                        sst["annual_sst_C"],
                        outdir / f"earth_reference_{int(cells)}cells_sst.png",
                        title="NOAA OISST v2 annual sea-surface temperature (C)",
                        cmap=render.TEMPERATURE_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-2.0,
                        vmax=32.0,
                    )
                    sst_seasons_png = _render_seasonal_panels(
                        grid,
                        sst["seasonal_sst_C"],
                        outdir / f"earth_reference_{int(cells)}cells_sst_seasons.png",
                        title="NOAA OISST sea-surface temperature (C)",
                        cmap=render.TEMPERATURE_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-2.0,
                        vmax=32.0,
                    )
                    ice_png = _render_scalar_map(
                        grid,
                        sst["annual_sea_ice_concentration_pct"],
                        outdir / f"earth_reference_{int(cells)}cells_sea_ice.png",
                        title="NOAA OISST annual sea-ice concentration (%)",
                        cmap="Blues",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=100.0,
                    )
                    for label, path in (
                        ("Annual SST", sst_png),
                        ("Seasonal SST", sst_seasons_png),
                        ("Sea Ice", ice_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))
            else:
                warnings.append("NOAA OISST files missing; SST reference skipped")

        if config.include_koppen:
            if koppen_available():
                kg = sample_koppen_to_grid(grid, elevation_m=elevation)
                koppen_class = np.asarray(kg["koppen_class"], dtype=np.int16)
                koppen_major = np.asarray(kg["koppen_major_class"], dtype=np.int16)
                biome_proxy = np.asarray(kg["biome_class_proxy"], dtype=np.int16)
                koppen_source_id = str(kg["koppen_source_id"])
                extra_arrays.update({
                    "earth.koppen_class": koppen_class,
                    "earth.koppen_major_class": koppen_major,
                    "earth.biome_class_proxy": biome_proxy,
                })
                metrics["koppen_class"] = categorical_fraction_metrics(
                    grid,
                    koppen_class,
                    prefix="koppen_class",
                )
                metrics["koppen_major"] = categorical_fraction_metrics(
                    grid,
                    koppen_major,
                    prefix="koppen_major",
                )
                metrics["biome_proxy"] = categorical_fraction_metrics(
                    grid,
                    biome_proxy,
                    prefix="biome_proxy",
                )
                metadata["source_ids"].append(koppen_source_id)
                metadata["koppen_class_names"] = kg["koppen_class_names"]
                metadata["koppen_native_resolution"] = kg["koppen_native_resolution"]
                metadata["koppen_period"] = kg["koppen_period"]
                if config.render_assets:
                    koppen_fine_png = _render_koppen_class_map(
                        grid,
                        koppen_class,
                        outdir / f"earth_reference_{int(cells)}cells_koppen_fine.png",
                        title=(
                            "GloH2O Koppen-Geiger 1991-2020 classes"
                            if koppen_source_id == "GLOH2O_KOPPEN_GEIGER"
                            else "Koppen-Geiger classes"
                        ),
                        width=int(config.width),
                        height=int(config.height),
                    )
                    koppen_png = _render_koppen_major_map(
                        grid,
                        koppen_major,
                        outdir / f"earth_reference_{int(cells)}cells_koppen_major.png",
                        title="Koppen-Geiger major classes",
                        width=int(config.width),
                        height=int(config.height),
                    )
                    biome_png = _render_biome_proxy_map(
                        grid,
                        biome_proxy,
                        outdir / f"earth_reference_{int(cells)}cells_biomes_from_koppen_proxy.png",
                        title="Biome proxy from Koppen-Geiger classes",
                        width=int(config.width),
                        height=int(config.height),
                    )
                    for label, path in (
                        ("Koppen Fine", koppen_fine_png),
                        ("Koppen Major", koppen_png),
                        ("Biome Proxy", biome_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))
            else:
                warnings.append("Koppen ASCII zip missing; climate-class reference skipped")

            if resolve_ecoregions_available():
                resolve_fields = sample_resolve_ecoregions_to_grid(
                    grid,
                    land_mask=elevation >= 0.0,
                )
                resolve_biome = np.asarray(
                    resolve_fields["resolve_biome_class"],
                    dtype=np.int16,
                )
                resolve_ecoregion = np.asarray(
                    resolve_fields["resolve_ecoregion_id"],
                    dtype=np.int32,
                )
                extra_arrays.update({
                    "earth.resolve_biome_class": resolve_biome,
                    "earth.resolve_ecoregion_id": resolve_ecoregion,
                })
                metrics["resolve_ecoregions_2017"] = resolve_ecoregion_metrics(
                    grid,
                    resolve_fields,
                    land_mask=elevation >= 0.0,
                )
                metadata["source_ids"].append("RESOLVE_ECOREGIONS_2017")
                metadata["resolve_biome_names"] = resolve_fields["resolve_biome_names"]
                metadata["resolve_ecoregion_names"] = resolve_fields[
                    "resolve_ecoregion_names"
                ]
                if config.render_assets:
                    resolve_png = _render_resolve_biome_map(
                        grid,
                        resolve_biome,
                        outdir / f"earth_reference_{int(cells)}cells_resolve_biomes.png",
                        title="RESOLVE Ecoregions 2017 terrestrial biomes",
                        width=int(config.width),
                        height=int(config.height),
                    )
                    asset_paths.append(str(resolve_png))
                    panel_paths.append(("RESOLVE Biomes", resolve_png))
            else:
                warnings.append(
                    "RESOLVE Ecoregions zip missing; true biome reference skipped"
                )

        if config.include_land_cover:
            try:
                land_cover = sample_esa_cci_land_cover_to_grid(
                    grid,
                    allow_download=bool(config.download_land_cover),
                )
            except FileNotFoundError as exc:
                warnings.append(str(exc))
            except Exception as exc:
                warnings.append(f"ESA CCI land-cover reference skipped: {exc}")
            else:
                extra_arrays.update({
                    "earth.esa_cci_lccs_class": land_cover["esa_cci_lccs_class"],
                    "earth.esa_cci_land_cover_broad_class": land_cover[
                        "esa_cci_land_cover_broad_class"
                    ],
                })
                metrics["esa_cci_land_cover"] = esa_cci_land_cover_metrics(
                    grid,
                    land_cover,
                    elevation_m=elevation,
                )
                metadata["source_ids"].append("ESA_CCI_LAND_COVER")
                metadata["esa_cci_land_cover"] = land_cover.get("metadata", {})
                metadata["land_cover_broad_names"] = LAND_COVER_BROAD_NAMES
                if config.render_assets:
                    lc_png = _render_scalar_map(
                        grid,
                        land_cover["esa_cci_land_cover_broad_class"],
                        outdir / f"earth_reference_{int(cells)}cells_land_cover_broad.png",
                        title="ESA CCI/C3S land cover broad classes",
                        cmap=LAND_COVER_BROAD_CMAP,
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=float(len(LAND_COVER_BROAD_COLORS) - 1),
                    )
                    asset_paths.append(str(lc_png))
                    panel_paths.append(("Land Cover", lc_png))

        if config.include_noaa_psl:
            if noaa_psl_available():
                ncep = sample_noaa_psl_to_grid(
                    grid,
                    start_year=int(config.climatology_start_year),
                    end_year=int(config.climatology_end_year),
                )
                extra_arrays.update({
                    "earth.monthly_wind_u10_v10": ncep["monthly_wind_u10_v10"],
                    "earth.seasonal_wind_u10_v10": ncep["seasonal_wind_u10_v10"],
                    "earth.monthly_slp_hPa": ncep["monthly_slp_hPa"],
                    "earth.seasonal_slp_hPa": ncep["seasonal_slp_hPa"],
                    "earth.annual_slp_hPa": ncep["annual_slp_hPa"],
                    "earth.seasonal_slp_anomaly_hPa": ncep[
                        "seasonal_slp_anomaly_hPa"
                    ],
                })
                metrics["noaa_psl_ncep"] = noaa_psl_metrics(ncep)
                metadata["source_ids"].append("NOAA_PSL_NCEP_NCAR_REANALYSIS_1")
                metadata["noaa_psl_climatology_years"] = [
                    int(config.climatology_start_year),
                    int(config.climatology_end_year),
                ]
                if config.render_assets:
                    wind_speed = np.linalg.norm(ncep["seasonal_wind_u10_v10"], axis=2)
                    wind_png = _render_seasonal_panels(
                        grid,
                        wind_speed,
                        outdir / f"earth_reference_{int(cells)}cells_wind_seasons.png",
                        title="NOAA PSL 10m wind speed (m/s)",
                        cmap="viridis",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=12.0,
                    )
                    pressure_png = _render_seasonal_panels(
                        grid,
                        ncep["seasonal_slp_anomaly_hPa"],
                        outdir / f"earth_reference_{int(cells)}cells_seasonal_pressure_anomaly.png",
                        title="NOAA PSL sea-level pressure anomaly (hPa)",
                        cmap="coolwarm",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=-12.0,
                        vmax=12.0,
                    )
                    for label, path in (
                        ("Seasonal Wind", wind_png),
                        ("Seasonal SLP Anomaly", pressure_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))
            else:
                warnings.append(
                    "NOAA PSL NCEP wind/pressure files missing; wind reference skipped"
                )

        if config.include_ocean_currents:
            primary_current_written = False
            try:
                oscar = sample_oscar_monthly_current_to_grid(
                    grid,
                    allow_download=bool(config.download_oscar),
                    level=OSCAR_ARCGIS_LEVEL,
                )
            except FileNotFoundError as exc:
                warnings.append(str(exc))
            except Exception as exc:
                warnings.append(f"NASA/JPL OSCAR monthly current reference skipped: {exc}")
            else:
                annual_uv = np.asarray(oscar["annual_surface_current_u_v"], dtype=np.float32)
                direction = (
                    np.degrees(np.arctan2(annual_uv[:, 0], annual_uv[:, 1])) % 360.0
                ).astype(np.float32)
                extra_arrays.update({
                    "earth.surface_current_u_v": annual_uv,
                    "earth.monthly_surface_current_u_v": oscar[
                        "monthly_surface_current_u_v"
                    ],
                    "earth.seasonal_surface_current_u_v": oscar[
                        "seasonal_surface_current_u_v"
                    ],
                    "earth.annual_surface_current_u_v": annual_uv,
                    "earth.monthly_surface_current_speed_m_s": oscar[
                        "monthly_surface_current_speed_m_s"
                    ],
                    "earth.seasonal_surface_current_speed_m_s": oscar[
                        "seasonal_surface_current_speed_m_s"
                    ],
                    "earth.annual_surface_current_speed_m_s": oscar[
                        "annual_surface_current_speed_m_s"
                    ],
                    "earth.annual_surface_current_direction_deg": direction,
                })
                metrics["nasa_jpl_oscar_monthly"] = oscar_monthly_current_metrics(
                    grid, oscar)
                metadata["source_ids"].append("NASA_JPL_OSCAR")
                metadata["oscar_arcgis"] = oscar.get("metadata", {})
                primary_current_written = True
                if config.render_assets:
                    oscar_png = _render_scalar_map(
                        grid,
                        oscar["annual_surface_current_speed_m_s"],
                        outdir / f"earth_reference_{int(cells)}cells_current_speed.png",
                        title="NASA/JPL OSCAR annual surface current speed (m/s)",
                        cmap="magma",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=1.0,
                    )
                    oscar_seasonal_png = _render_seasonal_panels(
                        grid,
                        oscar["seasonal_surface_current_speed_m_s"],
                        outdir / f"earth_reference_{int(cells)}cells_current_speed_seasons.png",
                        title="NASA/JPL OSCAR seasonal surface current speed (m/s)",
                        cmap="magma",
                        width=int(config.width),
                        height=int(config.height),
                        vmin=0.0,
                        vmax=1.0,
                    )
                    for label, path in (
                        ("Current Speed", oscar_png),
                        ("Seasonal Current Speed", oscar_seasonal_png),
                    ):
                        asset_paths.append(str(path))
                        panel_paths.append((label, path))

            if noaa_aoml_drifter_current_available():
                try:
                    currents = sample_noaa_aoml_drifter_current_to_grid(grid)
                except Exception as exc:
                    warnings.append(
                        f"NOAA/AOML drifter current reference skipped: {exc}"
                    )
                else:
                    aoml_arrays = {
                        "earth.aoml_surface_current_u_v": currents["surface_current_u_v"],
                        "earth.aoml_annual_surface_current_speed_m_s": currents[
                            "annual_surface_current_speed_m_s"
                        ],
                        "earth.aoml_annual_surface_current_direction_deg": currents[
                            "annual_surface_current_direction_deg"
                        ],
                    }
                    if not primary_current_written:
                        aoml_arrays.update({
                            "earth.surface_current_u_v": currents["surface_current_u_v"],
                            "earth.annual_surface_current_speed_m_s": currents[
                                "annual_surface_current_speed_m_s"
                            ],
                            "earth.annual_surface_current_direction_deg": currents[
                                "annual_surface_current_direction_deg"
                            ],
                        })
                    extra_arrays.update(aoml_arrays)
                    metrics["noaa_aoml_drifter_current_v3"] = (
                        noaa_aoml_drifter_current_metrics(grid, currents)
                    )
                    metadata["source_ids"].append(
                        "NOAA_AOML_DRIFTER_CURRENT_CLIMATOLOGY_V3"
                    )
                    if config.render_assets:
                        aoml_name = (
                            f"earth_reference_{int(cells)}cells_aoml_current_speed.png"
                            if primary_current_written
                            else f"earth_reference_{int(cells)}cells_current_speed.png"
                        )
                        current_png = _render_scalar_map(
                            grid,
                            currents["annual_surface_current_speed_m_s"],
                            outdir / aoml_name,
                            title="NOAA/AOML annual surface current speed (m/s)",
                            cmap="magma",
                            width=int(config.width),
                            height=int(config.height),
                            vmin=0.0,
                            vmax=1.0,
                        )
                        asset_paths.append(str(current_png))
                        panel_paths.append((
                            "AOML Current Speed" if primary_current_written else "Current Speed",
                            current_png,
                        ))
            else:
                warnings.append("lerc package missing; ocean-current reference skipped")

        if config.render_assets:
            contact = _render_contact_sheet(
                panel_paths,
                outdir / f"earth_reference_{int(cells)}cells_contact_sheet.png",
            )
            if contact is not None:
                asset_paths.append(str(contact))

        arrays_path = _write_arrays(
            outdir / f"earth_reference_{int(cells)}cells.npz",
            grid,
            elevation,
            metrics,
            extra_arrays=extra_arrays,
            extra_metadata=metadata,
        )
        entries.append({
            "cells": int(cells),
            "source_ids": metadata["source_ids"],
            "arrays": str(arrays_path),
            "assets": asset_paths,
            "metrics": metrics,
        })

    manifest = build_source_manifest()
    manifest_path.write_text(json.dumps(manifest, indent=2, default=_json_default))
    if config.manifest_out is not None:
        write_source_manifest(Path(config.manifest_out))

    completed_sources = {
        source_id
        for entry in entries
        for source_id in entry.get("source_ids", [])
    }
    next_reference_sources = [
        source
        for source in (
            "NASA_JPL_OSCAR",
            "ESA_CCI_LAND_COVER",
            "MODIS_MCD12C1",
        )
        if source not in completed_sources
    ]
    summary = {
        "schema": SCHEMA,
        "source_manifest": str(manifest_path),
        "data_manifest": str(config.manifest_out) if config.manifest_out else None,
        "cells": [int(c) for c in config.cells],
        "entry_count": len(entries),
        "entries": entries,
        "warnings": warnings,
        "next_reference_sources": next_reference_sources,
    }
    (outdir / "earth_climate_reference_summary.json").write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    return summary
