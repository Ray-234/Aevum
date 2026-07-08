# Earth Climate Reference Calibration Plan

Status: active planning document; R0/R1/R2/R3/R4 Earth references active
Owner: climate / biome calibration
Last updated: 2026-07-06

This plan adds a real-Earth reference track for calibrating Aevum's generated
climate and biome maps.  It is intentionally separate from plate/terrain
generation: the accepted plate system remains frozen while climate outputs are
compared against observed Earth fields.

## Goal

Build a repeatable reference library that can:

- Download or ingest public Earth datasets.
- Resample them onto Aevum's spherical cell grid and standard raster outputs.
- Produce side-by-side maps with generated worlds.
- Compute calibration metrics for terrain, winds, currents, temperature,
  precipitation, climate classes, and biomes.

The purpose is not to force generated worlds to look exactly like Earth.  The
purpose is to keep first-order Earth physics visible enough that an Earthlike
layout produces Earthlike climate logic.

## Candidate Reference Datasets

Primary data sources:

- Topography and bathymetry:
  - NOAA/NCEI ETOPO Global Relief Model / ETOPO 2022:
    https://www.ncei.noaa.gov/products/etopo-global-relief-model
    https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem%3Aetopo_2022
- Land monthly temperature and precipitation:
  - WorldClim v2.1 historical climate, 1970-2000 normals:
    https://www.worldclim.org/data/worldclim21.html
- Global monthly atmospheric fields:
  - ERA5 monthly averaged single levels, 1940-present:
    https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means
  - NOAA PSL NCEP/NCAR Reanalysis 1 as a no-account fallback:
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html
- Surface ocean currents:
  - NASA/JPL PO.DAAC OSCAR surface currents v2.0:
    https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0
    https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_1deg
- Koppen-Geiger climate classes:
  - Beck et al. / GloH2O historical and future Koppen-Geiger maps:
    https://www.gloh2o.org/koppen/
    https://www.nature.com/articles/sdata2018214
- Biome / land-cover references:
  - ESA CCI/Copernicus Land Cover:
    https://esa-landcover-cci.org/
    https://climate.esa.int/en/projects/land-cover/
  - MODIS MCD12C1 global land-cover type:
    https://modis.gsfc.nasa.gov/data/dataprod/mod12.php
    https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MCD12C1
  - RESOLVE / Ecoregions 2017 biome polygons:
    https://ecoregions.appspot.com/

## Data Policy

- Prefer low or moderate resolution products first.  Calibration only needs
  global structure at Aevum cell scales, not native 300 m or 1 km resolution.
- Cache raw downloads under `data/reference/earth_climate/`.
- Store derived grid-cell arrays as compressed `.npz` with provenance metadata.
- Keep any source requiring an account, license click-through, or API token
  optional.  ERA5 is high-quality but should not block the pipeline; NOAA PSL
  and WorldClim can provide fallback climatology.
- Do not commit large raw datasets by default.  Commit small manifests,
  checksums, scripts, and derived low-resolution fixtures only.

## Reference Pipeline

### R0 - Source Ledger

Deliverables:

- `data/reference/earth_climate/source_manifest.json`
- One row per dataset: source URL, license/citation note, variables, temporal
  baseline, native resolution, local cache path, checksum if downloaded.

Acceptance:

- Every reference map can be traced back to a public dataset and version.

### R1 - Earth Grid Adapter

Deliverables:

- A sampler that maps lat/lon rasters or NetCDF grids to `SphereGrid` cells.
- Standard derived arrays:
  - `earth.elevation_m`
  - `earth.land_mask`
  - `earth.monthly_temperature_C`
  - `earth.monthly_precip_mm`
  - `earth.seasonal_temperature_C`
  - `earth.seasonal_precip_mm_yr_equiv`
  - `earth.wind_u10_v10`
  - `earth.surface_current_u_v`
  - `earth.koppen_class`
  - `earth.biome_class`

Acceptance:

- Cell arrays have shape `(n_cells,)`, `(4, n_cells)`, or `(12, n_cells)` as
  appropriate and are finite except for documented ocean/land masks.
- Antimeridian continuity is preserved.

### R2 - Reference Rendering

Deliverables:

- `earth_elevation.png`
- `earth_currents.png`
- `earth_wind_seasons.png`
- `earth_temperature_seasons.png`
- `earth_precip_seasons.png`
- `earth_koppen.png`
- `earth_biomes.png`

Acceptance:

- Earth maps use the same projection, extents, and comparable color scales as
  generated Aevum maps.

### R3 - Calibration Metrics

Terrain metrics:

- Hypsometry: land mean/p50/p90, ocean mean/p90 depth.
- Shelf/deep-ocean fraction and nearshore depth envelope.
- Mountain-belt continuity proxy from elevation gradients.

Climate metrics:

- Zonal mean temperature curve and adjacent-latitude jump.
- Land/ocean seasonal temperature amplitude contrast.
- Seasonal ITCZ/rain-belt migration proxy.
- Subtropical dry belt latitude and continuity.
- Mid-latitude storm-track wet-coast signal.
- Monsoon index over large heated continents.
- Rain-shadow/coastal-wetness asymmetry.

Ocean-current metrics:

- Current speed distribution.
- Western-boundary vs eastern-boundary thermal/asymmetry proxy.
- Upwelling/cold-current dry-coast co-location.

Biome metrics:

- Class fractions by latitude band.
- Desert/grassland/forest/tropical/tundra/ice envelopes.
- Agreement between generated biome thresholds and Koppen/land-cover reference
  at coarse class level.

Acceptance:

- Earth reference metrics are produced first as baseline values.
- Generated Earthlike worlds are compared to Earth envelopes, not one exact
  pixel arrangement.

### R4 - Calibration Loops

Use the reference metrics to tune climate parameters in this order:

1. Temperature: annual mean, lapse-rate expression, seasonal amplitude.
2. Moisture access: coastal wetness, interior drying, rain-shadow strength.
3. Seasonal precipitation: ITCZ migration, monsoon potential, dry-season length.
4. Ocean influence: warm/cold-current thermal and precipitation asymmetry.
5. Biomes: class thresholds, dry/cold seasonal stress, tropical/temperate split.

Each tuning pass should generate:

- One Earth reference contact sheet.
- One generated six-world contact sheet.
- Metrics JSON with before/after deltas.

## Current State

R0 source manifest and executable R1/R2/R3 adapters now exist in:

- `aevum.diagnostics.earth_climate_reference`
- CLI command: `aevum earth-climate-reference`
- Persistent manifest target:
  `data/reference/earth_climate/source_manifest.json`
- Current generated reference assets:
  `out_earth_climate_reference_r4_20260705/`
- Current OSCAR monthly-current reference assets:
  `out_earth_climate_reference_r5_oscar_20260706/`
- Current ESA CCI land-cover reference assets:
  `out_earth_climate_reference_r6_landcover_20260706/`
- Current ETOPO 2022 cross-check assets:
  `out_earth_climate_reference_r4_etopo2022_render_20260706/`
  `out_earth_climate_reference_r4_etopo2022_24000_render_20260706/`

Implemented source adapters:

- ETOPO5 local relief: `earth.elevation_m`, `earth.land_mask`, hypsometry.
- ETOPO 2022 lightweight topography/bathymetry cross-check:
  `earth.etopo2022_elevation_m`, `earth.etopo2022_land_mask`, and
  `earth.etopo2022_minus_etopo5_m`.  The cached file is the NOAA THREDDS
  60 arc-second bedrock relief product sampled every 30 cells through OPeNDAP,
  yielding a 0.5 degree derived cache under
  `data/reference/earth_climate/raw/etopo_2022/`.
- WorldClim v2.1 10 arc-minute monthly land climate:
  `earth.monthly_temperature_C`, `earth.monthly_precip_mm`,
  `earth.seasonal_temperature_C`, `earth.seasonal_precip_mm_yr_equiv`,
  `earth.annual_temperature_C`, `earth.annual_precip_mm`, and
  `earth.dry_month_count`.
- GloH2O / Beck et al. Koppen-Geiger V3, default 1991-2020 0.1 degree
  GeoTIFF:
  `earth.koppen_class`, `earth.koppen_major_class`, and
  `earth.biome_class_proxy`.
- Kottek/Rubel 0.5 degree Koppen-Geiger ASCII remains available as fallback
  when the GloH2O zip is absent.
- RESOLVE Ecoregions 2017:
  `earth.resolve_biome_class` and `earth.resolve_ecoregion_id`.
- NOAA OISST v2 monthly long-term mean, 1991-2020:
  `earth.monthly_sst_C`, `earth.seasonal_sst_C`, `earth.annual_sst_C`,
  `earth.monthly_sea_ice_concentration_pct`,
  `earth.seasonal_sea_ice_concentration_pct`, and
  `earth.annual_sea_ice_concentration_pct`.
- NOAA PSL NCEP/NCAR Reanalysis 1 monthly climatology, default 1991-2020:
  `earth.monthly_wind_u10_v10`, `earth.seasonal_wind_u10_v10`,
  `earth.monthly_slp_hPa`, `earth.seasonal_slp_hPa`,
  `earth.annual_slp_hPa`, and `earth.seasonal_slp_anomaly_hPa`.
- ESA CCI / C3S 2020 land cover through Planetary Computer STAC rendered
  categorical previews:
  `earth.esa_cci_lccs_class` and
  `earth.esa_cci_land_cover_broad_class`.  This is a coarse cross-check
  against land-cover envelopes, not a replacement for full 300 m land-use
  analysis.
- NASA/JPL PO.DAAC OSCAR monthly surface-current climatology through the public
  ArcGIS Vector-UV LERC2D cache, sampled from 240 monthly slices in 2001-2020:
  `earth.monthly_surface_current_u_v`,
  `earth.seasonal_surface_current_u_v`,
  `earth.annual_surface_current_u_v`,
  `earth.monthly_surface_current_speed_m_s`,
  `earth.seasonal_surface_current_speed_m_s`, and
  `earth.annual_surface_current_speed_m_s`.  The compact level-1 monthly
  climatology cache is under `data/reference/earth_climate/raw/oscar/`.
- NOAA/AOML annual near-surface drifter current climatology v3 remains as a
  no-account cross-check:
  `earth.aoml_surface_current_u_v`,
  `earth.aoml_annual_surface_current_speed_m_s`, and
  `earth.aoml_annual_surface_current_direction_deg`.  It is still promoted to
  the generic `earth.surface_current_u_v` fields only when OSCAR is absent.

Generated R4/R5 maps at 8000 and 24000 Aevum-cell sampling density:

- `earth_reference_*cells_elevation.png`
- `earth_reference_*cells_etopo2022_elevation.png`
- `earth_reference_*cells_etopo2022_minus_etopo5.png`
- `earth_reference_*cells_temperature.png`
- `earth_reference_*cells_precip.png`
- `earth_reference_*cells_temperature_seasons.png`
- `earth_reference_*cells_precip_seasons.png`
- `earth_reference_*cells_sst.png`
- `earth_reference_*cells_sst_seasons.png`
- `earth_reference_*cells_sea_ice.png`
- `earth_reference_*cells_koppen_fine.png`
- `earth_reference_*cells_koppen_major.png`
- `earth_reference_*cells_biomes_from_koppen_proxy.png`
- `earth_reference_*cells_resolve_biomes.png`
- `earth_reference_*cells_wind_seasons.png`
- `earth_reference_*cells_seasonal_pressure_anomaly.png`
- `earth_reference_*cells_current_speed.png`
- `earth_reference_*cells_current_speed_seasons.png`
- `earth_reference_*cells_aoml_current_speed.png`
- `earth_reference_*cells_land_cover_broad.png`
- `earth_reference_*cells_contact_sheet.png`

Current 8000/24000 reference metrics:

- Land fraction: `0.287` / `0.286`
- Mean land elevation: `810 m` / `820 m`
- Mean ocean depth: `3680 m` / `3668 m`
- WorldClim land annual temperature mean: `9.45 C` / `9.45 C`
- WorldClim land annual precipitation mean: `767 mm/yr` / `761 mm/yr`
- WorldClim median land dry-month count: `6` / `6`
- NOAA PSL seasonal wind speed p90: `7.0 m/s` / `7.0 m/s`
- NOAA PSL seasonal sea-level-pressure anomaly absolute p90: `4.1 hPa` /
  `4.1 hPa`
- NOAA OISST annual ocean SST mean: `18.3 C` / `18.3 C`
- NOAA OISST tropical ocean SST mean: `27.0 C` / `27.0 C`
- NOAA OISST seasonal SST amplitude p90: `5.45 C` / `5.48 C`
- NOAA OISST annual sea-ice area fraction above 15% concentration:
  `0.089` / `0.089` of valid ocean area.
- RESOLVE assigned land-area fraction: `0.988` / `0.986`
- RESOLVE observed biome count: `14` / `14`
- 24000-cell sampled RESOLVE ecoregion ids observed: `616`
- 24000-cell sampled GloH2O Koppen classes observed: `27`
- NOAA/AOML drifter current valid area fraction: `0.628` / `0.627`
- NOAA/AOML drifter current speed p90: `0.27 m/s` / `0.26 m/s`
- NASA/JPL OSCAR current valid area fraction: `0.693` / `0.695`
- NASA/JPL OSCAR annual current speed p50: `0.051 m/s` / `0.051 m/s`
- NASA/JPL OSCAR annual current speed p90: `0.192 m/s` / `0.190 m/s`
- NASA/JPL OSCAR monthly current speed p90: `0.223 m/s` / `0.225 m/s`
- NASA/JPL OSCAR seasonal current-speed amplitude p90:
  `0.163 m/s` / `0.162 m/s`
- ESA CCI land-cover valid area fraction: `1.000` / `1.000`
- ESA CCI observed LCCS class count: `26` / `26`
- ESA CCI broad water fraction: `0.718` / `0.718`
- ESA CCI broad forest fraction: `0.084` / `0.085`
- ESA CCI broad cropland fraction: `0.041` / `0.042`
- ESA CCI broad grass/shrub fraction: `0.068` / `0.065`
- ESA CCI land-cover water-vs-ETOPO5-ocean mismatch area:
  `0.013` / `0.013`

Current ETOPO 2022 8000/24000-cell cross-check metrics:

- ETOPO 2022 land fraction: `0.276` / `0.277`, versus ETOPO5
  `0.287` / `0.286` (`-0.0116` / `-0.0087` delta).
- ETOPO 2022 mean land elevation: `644 m` / `644 m`, versus ETOPO5
  `810 m` / `820 m` (`-166 m` / `-176 m` delta).
- ETOPO 2022 mean ocean depth: `3633 m` / `3634 m`, versus ETOPO5
  `3680 m` / `3668 m` (`-46 m` / `-34 m` delta).
- Mean absolute ETOPO2022-minus-ETOPO5 elevation delta:
  `247 m` / `251 m`; p95 absolute delta: `976 m` / `988 m`.
- ETOPO5/ETOPO2022 land-ocean classification mismatch area:
  `0.026` / `0.027`.

Remaining initial dataset set:

1. MODIS MCD12C1 remains an optional secondary land-cover cross-check.  ESA CCI
   now covers the primary land-cover reference requirement.
2. Decide whether to promote ETOPO 2022 from an uncertainty cross-check to the
   primary terrain reference after the climate gates are compared against both
   ETOPO5 and ETOPO 2022 envelopes.

Next executable step:

- Treat the Earth reference pipeline as sufficient for the next climate fitting
  pass.  MODIS MCD12C1 can be added later as a secondary land-cover audit if
  ESA CCI and RESOLVE disagree in a way that blocks biome tuning.
- Keep OSCAR as the current primary current-speed calibration source for new
  R5 Earth references.  AOML remains useful for checking that the OSCAR
  level-1 climatology is not under- or over-scaling the annual speed envelope.
- R6 has now been used to compare the current rendered C4j six-world terminal
  climate/biome replay against Earth envelopes:
  `out_earth_climate_comparison_c4j1_render_r6_20260706/` reports
  `earthlike flagged: 0`, and
  `out_earth_climate_fitting_c4j1_render_r6_20260706/` reports guardrail
  verdict `pass` with `0` failures and `0` warnings.
- Use the ETOPO 2022 cross-check metrics as an uncertainty band for Earth
  hypsometry and bathymetry gates instead of silently treating ETOPO5 as exact.
