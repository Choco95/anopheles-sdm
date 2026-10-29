# Data guide

## Nine-species SDM

The following filenames and columns are required by the retained code. Use GeoTIFF extension `.tif` or `.tiff`. The code checks raster stems case-insensitively, but portable projects should preserve the documented spelling.

### Occurrences

Place these files in the configured `OCC_DIR` (recommended `data/raw/occurrences/`).

| Species | CSV filename |
| --- | --- |
| funestus | `funestus.csv` |
| arabiensis | `Arabiensis.csv` |
| coluzzii | `Coluzzi.csv` |
| coustani | `Coustani.csv` |
| merus | `Merus.csv` |
| quadriannulatus | `quadra.csv` |
| stephensi | `Stephensi.csv` |
| moucheti | `moucheti.csv` |
| nili | `nili.csv` |

Each file must contain exact columns `Latitude`, `Longitude`, `Species`. Coordinates are expected in WGS84 longitude/latitude. Additional source fields may be kept in the original file, but the script's cleaning function retains only these three core columns plus its audit fields. Keep the original files to retain dates, sources and sampling attributes. Exact accepted species labels are in `ACCEPTED_LABELS`; file names do not validate the taxonomic identity of rows.

### Environmental inputs

`ENV_DIR/bio_clim_var/` must contain exactly nineteen identifiable, single-band BIO1–BIO19 rasters. Filenames must contain a `bioN`, `bio_N` or `bio-N` identifier accepted by the script. Do not put duplicate product versions in that directory.

The following additional rasters belong directly in `ENV_DIR`:

| Standardised SDM name | Required source filename | Data kind |
| --- | --- | --- |
| lulc | `ESA_LULC_Africa.tif` | Categorical |
| human_footprint | `Human_Footprint.tif` | Continuous |
| humidity | `Humidity.tif` | Continuous |
| ndvi | `NDVI.tif` | Continuous |
| population_density | `populationDensity.tif` | Continuous |
| slope | `Slope.tif` | Continuous |
| solar_radiation | `Solar Radiation.tif` | Continuous |
| elevation | `wc2.1_5m_elev.tif` | Continuous |
| wind_speed | `WindSpeed.tif` | Continuous |
| aspect | `Aspect.tif` | Angle, converted to northness/eastness |
| distance_to_water | `DistanceToWater.tif` | Continuous |
| twi | `TWI.tif` | Continuous |
| water | `WaterRasterized.tif` | Binary |

Source product, version, units and time coverage are not established by this table. Fill `metadata/covariates.csv` from source documentation and raster metadata. Native resolution and the aligned working resolution are different metadata fields. The code does not automatically average multiband inputs into temporal summaries.

## Separate five-species overlap

The overlap script reconstructs fourteen predictors: `bio2`, `bio6`, `bio14`, `bio15`, `Humidity`, `Elevation`, `Slope`, `AspectEast`, `AspectNorth`, `TWI`, `DistanceToWater`, `LULC`, `PopulationDensity`, `HumanFootprint`.

Its source-file lookup differs from the nine-species SDM:

| Overlap filename | Nine-species filename | What to verify |
| --- | --- | --- |
| `Elevation.tif` | `wc2.1_5m_elev.tif` | Same actual elevation dataset if linked |
| `LULC.tif` | `ESA_LULC_Africa.tif` | Same category coding/product |
| `HumanFootprint.tif` | `Human_Footprint.tif` | Same product/units |
| `PopulationDensity.tif` | `populationDensity.tif` | Lookup ignores case; same product still required |
| `AspectNorth.tif` / `AspectEast.tif` | Computed northness/eastness outputs | Confirm compatible transformation, flat-cell treatment and grid |

Other non-BIO files are `Humidity.tif`, `Slope.tif`, `TWI.tif`, `DistanceToWater.tif`. If both aspect-component rasters are unavailable, the script reads `Aspect.tif`, aligns it and then computes components; that differs from the nine-species script's transform-before-alignment sequence.

The `inputs` table must supply occurrence CSV, continuous map and binary map paths for arabiensis, coluzzii, coustani, funestus and gambiae. The retained overlap paths expect existing `ensemble_map.tif` and `binary_map.tif`; the nine-species script writes `ensemble_suitability.tif` and `binary_suitability.tif`. Edit paths deliberately rather than assuming these names refer to equivalent outputs. Record each map's actual source run in `metadata/overlap_inputs.csv`.

All maps must be single-band, with identical grid, CRS and extent. Binary maps contain only 0/1 plus NA. The overlap code takes the common complete domain across predictors and maps. It does not derive the missing *An. gambiae* maps.

## Data access and checksums

Download locations, dataset versions and extraction dates must be completed from the actual source records. No occurrence data, rasters, credentials or fitted models are bundled.

For a local integrity manifest in R, use the actual configured input roots:

```r
files <- sort(list.files("data/raw", recursive = TRUE, full.names = TRUE))
files <- files[!file.info(files)$isdir & !grepl("README\\.md$", files)]
manifest <- data.frame(path = files, md5 = unname(tools::md5sum(files)))
write.csv(manifest, "logs/input_checksums.csv", row.names = FALSE)
```

If data are stored outside `data/raw`, enumerate those directories instead. The example uses MD5 for local change detection, not security verification. Archive this manifest with the run; it is excluded from routine Git commits under `logs/`.
