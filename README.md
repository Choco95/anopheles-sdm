# Anopheles species distribution modelling and niche overlap

R workflows and operating procedures for modelling *Anopheles* environmental suitability and comparing species distributions.

**Status: documentation and repository draft, 29 October 2026.** The two existing R scripts are included unchanged. Input datasets and fitted results are not included. No models were executed while preparing this repository, and this is not yet a verified reproduction of the reported study results.

Prepared for Marc Marin Agbodjan, icipe, in the context of IMPACTING WP1 reproducibility discussions. Project attribution, contributors and release details remain to be confirmed.

## Start here

- [SOP: full procedure and quality checks](docs/SOP.md)
- [Data requirements and exact filenames](docs/DATA_GUIDE.md)
- [Known differences and work remaining](docs/IMPLEMENTATION_STATUS.md)
- [GitHub setup](docs/GITHUB_SETUP.md)

## Included workflows

| Script | Species and purpose | Important boundary |
| --- | --- | --- |
| `scripts/Anopheles_Africa_SDM_updated.R` | Nine species; raster/occurrence audit, spatial validation, tuning, current suitability, optional future projection, geographic and PCA-env overlap | Starts in `RUN_MODE = "prepare"`; the nine-species list does not include *An. gambiae* |
| `scripts/Anopheles_niche_overlap_5_species.R` | *An. arabiensis*, *An. coluzzii*, *An. coustani*, *An. funestus*, *An. gambiae*; ten pairwise comparisons using existing maps | Does not fit SDMs; needs a separately supplied *An. gambiae* SDM and matching maps |

The nine-species script specifies *An. funestus*, *An. arabiensis*, *An. coluzzii*, *An. coustani*, *An. merus*, *An. quadriannulatus*, *An. stephensi*, *An. moucheti* and *An. nili*. Verify that these scripts match the version used for each figure before describing them as the final study implementation.

## Repository layout

| Location | Contents |
| --- | --- |
| `scripts/` | Two existing analysis scripts, package setup and a check that does not run analyses |
| `docs/` | SOP, data guide, implementation notes and GitHub instructions |
| `metadata/` | Covariate, occurrence and run-record templates; source checksums |
| `data/raw/occurrences/` | Local occurrence CSVs; data contents excluded from Git |
| `data/raw/environment/` | Local GeoTIFFs; data contents excluded from Git |
| `results/` | Local generated outputs; contents excluded from Git |
| `logs/` | Local run logs; contents excluded from Git |

## Software setup

Install R and RStudio, open `anopheles-sdm.Rproj`, then run:

```r
source("scripts/00_setup.R")
source("scripts/00_check_environment.R")
```

The setup script installs missing dependencies; the check script parses the analysis code without sourcing it and reports package availability. It does not check raster inputs or scientific validity. An exact R version and package lockfile have not been supplied.

For the first successful local execution, record the R version and create an environment lockfile using [renv](https://rstudio.github.io/renv/articles/renv.html). With this project open, use `renv::init(bare = TRUE)`, run `scripts/00_setup.R` in the activated project, then `renv::snapshot(type = "all")`. Commit the resulting lockfile and activation files. On later machines use `renv::restore()` instead of reinstalling current packages. Also record GDAL/GEOS/PROJ versions where relevant; an R lockfile does not install these external system libraries.

## Run the SDM audit first

The original scripts retain their Windows paths. In the top settings block of `Anopheles_Africa_SDM_updated.R`, replace the directory settings with:

```r
ROOT <- normalizePath(getwd(), winslash = "/")
ENV_DIR <- file.path(ROOT, "data", "raw", "environment")
OCC_DIR <- file.path(ROOT, "data", "raw", "occurrences")
OUT_DIR <- file.path(ROOT, "results", "sdm_run_001")
RUN_MODE <- "prepare"
```

Run from the project root in a fresh R session. Populate the input folders using [DATA_GUIDE.md](docs/DATA_GUIDE.md), save the settings, then run:

```r
source("scripts/Anopheles_Africa_SDM_updated.R")
```

Review `occurrence_summary.csv`, per-species exclusions, the raster audit and predictor coverage map. The preparation stage aligns rasters and audits records; it does not fit models.

After resolving the audit, change `RUN_MODE` to `"all"` and source the same script locally to fit the models. Use a new output directory for changed inputs, configuration, functions or package versions. The script can resume checkpoints for an unchanged run. Keep each invocation in a fresh R session because the fingerprint inspects functions in the global environment.

## Settings already present in the nine-species script

| Setting | Included value | Interpretation |
| --- | --- | --- |
| Algorithms | Maxnet, ranger RF, XGBoost, mgcv GAM | Implemented algorithms, not a claim about published analyses |
| Validation | 5 outer / 3 inner spatial folds | Review fold counts before accepting results |
| Spatial block size | 100 km | Provisional setting requiring assessment |
| Background | Up to 10,000 random available cells | Background cells are not observed absences |
| Calibration domain | Africa | Broad-domain baseline; justify or supply accessible-area polygons |
| Predictor screening | Correlation limit 0.70; VIF limit 5 | Selection occurs within training partitions |
| Ensemble weights | Training-derived TSS | Inspect saved weights and thresholds |
| Future configuration | Empty `FUTURE` table | No future projections run by default |

The script uses the geometry of its reference BIO raster. The phrase “~5 km” in its header does not verify an input's actual cell size. Record native resolution, final resolution, CRS and temporal coverage in `metadata/covariates.csv`.

## Run the separate five-species overlap analysis

In a fresh R session, edit `niche_out_dir`, `environment_root`, and every path in the `inputs` table of `scripts/Anopheles_niche_overlap_5_species.R`. The overlap script uses different environmental filenames from the SDM script; consult the data guide. Point it to compatible single-band continuous and binary maps for all five species, including *An. gambiae*.

```r
source("scripts/Anopheles_niche_overlap_5_species.R")
```

This produces `pairwise_niche_overlap.csv`, occurrence audits, settings and three metric heatmaps. Map export is off by default; set `save_overlap_maps <- TRUE` for the four-class maps.

The retained script contains an exploratory SVM cutoff of 0.80 and a two-of-three classification rule. Sørensen and centroid cutoffs are unset. Treat these as recorded script settings, not established biological thresholds. Report continuous metrics and do not interpret unresolved classifications as evidence of no overlap. See the SOP for definitions and limitations.

## Expected SDM outputs

| Output | Purpose |
| --- | --- |
| `occurrence_summary.csv` | Record retention across species |
| `environment_current/input_raster_audit.csv` | Source raster geometry audit |
| `<species>/spatial_CV_summary.csv` | Outer spatial validation summary |
| `<species>/current_maps/ensemble_suitability.tif` | Relative environmental suitability |
| `<species>/current_maps/binary_suitability.tif` | Binary map at the saved species threshold |
| `<species>/current_maps/algorithm_disagreement_SD.tif` | Between-algorithm disagreement diagnostic |
| `<species>/current_maps/environmental_novelty_flag.tif` | Univariate range novelty diagnostic |
| `overlap_current/pairwise_geographic_overlap_36_pairs.csv` | Nine-species geographic overlap |
| `niche_overlap/PCA_env_overlap_36_pairs.csv` | Separate PCA-env overlap analysis |
| `configuration.rds`, `sessionInfo.txt`, `input_file_manifest.csv` | Run configuration and environment record |

## Data, interpretation and release

Keep original inputs unchanged. Document source URLs, dataset versions, extraction dates, units, transformations and redistribution terms. Download instructions and data citations must be completed from the actual files; source product names cannot be inferred reliably from filenames alone.

Suitability is not mosquito abundance, observed co-occurrence, malaria incidence or transmission probability. The repository does not currently implement SEIR coupling.

No reuse license or publication DOI has been assigned here. Confirm ownership and the intended license before release. Dataset terms remain separate from any eventual code license. Fill in contributors and the manuscript citation when agreed.
