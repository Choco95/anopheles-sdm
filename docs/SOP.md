# Standard Operating Procedure

## Document control

| Field | Value |
| --- | --- |
| SOP identifier | ANOPHELES-SDM-001 |
| Title | Preparation, execution, quality assurance and archiving of Anopheles SDM and niche-overlap analyses |
| Version | 0.1 |
| Date | 29 October 2026 |
| Status | Draft for technical review; not an approved or validated protocol |
| Prepared for | Marc Marin Agbodjan, icipe |
| Technical reviewer / approval date | To be assigned / not yet approved |
| Source code | Two retained R scripts listed in `metadata/source_provenance.csv` |

## 1. Purpose and scope

Provide a traceable procedure for organising environmental and occurrence data, auditing inputs, running species distribution models, evaluating predictions, calculating overlap and retaining the information needed to reproduce a run.

This SOP describes the code currently included in the repository. It distinguishes the nine-species Africa SDM workflow from the later five-species overlap workflow. Their species lists, predictor naming, map filenames and overlap methods differ. The overlap analysis is not an automatic continuation of the nine-species script: it needs separately compatible inputs, including an *An. gambiae* model absent from that script.

Future projections are optional and require a populated scenario table. Transmission modelling and SEIR coupling fall outside this implementation.

## 2. Responsibilities

The analyst records data provenance, verifies inputs, configures the scripts, examines diagnostic outputs and preserves run records. A technical reviewer examines species definitions, calibration domain, validation design, predictor processing and interpretation before results are released. Assign actual people and dates in the run log; these roles do not imply institutional approval has already occurred.

## 3. Required materials

- R, RStudio and the packages declared in `scripts/00_setup.R`.
- The occurrence CSVs and environmental GeoTIFFs specified in `DATA_GUIDE.md`.
- A study boundary or the script's default Natural Earth Africa boundary; optional species calibration polygons when that mode is selected.
- For future runs, harmonised predictor sets for each GCM, SSP and period.
- For five-species overlap, five occurrence tables and ten corresponding continuous/binary maps from identified SDM runs.
- Local disk space sufficient for aligned Africa rasters, model checkpoints and outputs. Runtime and storage requirements must be measured locally.

## 4. Procedure

### 4.1 Establish the project and record the run

1. Open `anopheles-sdm.Rproj` and start a fresh R session.
2. Create a distinct run identifier and output directory. Enter operator, date, purpose, code commit and intended workflow in `metadata/run_log.csv`.
3. Record whether this is a preparation audit, full SDM run, future projection or five-species overlap run. A completed audit is not a completed modelling run.
4. Run the setup/check scripts. Record the exact R and package versions from the first successful local execution and create an `renv` lockfile as described in the README. Until then, package versions remain unpinned.

**Check:** the two analysis files parse and required packages load. Stop to resolve missing packages or syntax errors; do not source an analysis file simply to check its syntax, because each file executes its workflow when sourced.

### 4.2 Register and preserve source data

1. Place local source rasters and occurrences in the documented input folders, or explicitly configure their external locations.
2. Preserve originals. Do not overwrite source rasters with aligned versions or edit records in place.
3. Complete `metadata/covariates.csv` for each supplied raster. Record provider/product, version, URL, citation, download date, represented period, native CRS/resolution, final CRS/resolution, units, scale/offset, missing-value convention, preprocessing and license.
4. Complete `metadata/occurrence_sources.csv`: source/version, taxonomic scope, observation period, coordinate CRS, original and retained counts, and access terms.
5. Retain a checksum manifest of the actual inputs alongside the run archive. The SDM's input manifest records paths, size and modification time; its fingerprint does not hash every raster's bytes. It therefore does not replace content checksums.

**Check:** each input is linked to an identifiable source and transformation record. Unknown product versions, time windows or units remain unresolved items and must not be invented from a filename.

### 4.3 Configure the nine-species SDM

1. Set `ROOT`, `ENV_DIR`, `OCC_DIR`, `OUT_DIR` and `RUN_MODE` as shown in the README. Begin with `RUN_MODE <- "prepare"`.
2. Review `SPECIES` and `ACCEPTED_LABELS`. Retain the original CSV filenames and exact required column names. Document any taxonomic changes; complex or uncertain labels must not silently become sensu-stricto records.
3. Review all `CFG` settings. In the included version, the calibration domain is Africa, random background is requested up to 10,000 cells, spatial blocks are 100 km, and validation uses five outer and three inner folds. These are implementation settings to assess, not study-specific values proven by this SOP.
4. If using accessible-area polygons, provide each species polygon and set `calibration_mode = "polygons"`. Explain the domain choice in the run record.
5. Verify units, scale factors and aspect conventions before entering multipliers or offsets. Avoid applying a product scale factor twice.
6. Leave `FUTURE` empty for current-only analysis. When adding scenarios, use actual GCM names, SSPs, periods and input directories; replace all example identifiers.

**Check:** no placeholder paths or example scenario names are active. Commit configured code or preserve its exact copy and checksum with the run. Configuration changes require a new output directory.

### 4.4 Audit predictor geometry and coverage

1. Source the SDM script in preparation mode. The script loads nineteen BIO layers, reads additional predictors and builds an aligned environmental stack.
2. Inspect `environment_current/input_raster_audit.csv` and `aligned_manifest.csv`. Confirm each raster's CRS, extent, cell dimensions, band count and expected units.
3. Inspect `valid_predictor_coverage.png` and the valid mask. All-candidate completeness can reduce the prediction domain even when a predictor is later excluded.
4. Confirm categorical variables use nearest-neighbour alignment and continuous variables use the implemented continuous interpolation. Inspect coastal boundaries and missing-data edges.
5. In the nine-species script, aspect is transformed to northness/eastness before alignment and flat terrain is handled explicitly. Record that transformation.

**Check:** aligned geometry matches the reference raster and the final valid domain is scientifically acceptable. Do not describe the reference grid as exactly 5 km without reading its geometry. Record the measured grid instead.

### 4.5 Audit occurrence records

1. Review `occurrence_summary.csv` and, for every species, `occurrence_audit_all_rows.csv`, `excluded_occurrences_review.csv`, `original_label_counts.csv` and `clean_occurrences.csv`.
2. Examine exclusions for taxonomic labels, invalid or zero coordinates, records outside Africa, exact duplicates, missing predictors and repeated modelling cells.
3. Check counts and locations against the source data, especially near islands and coasts where the default coarse boundary can exclude valid points.
4. In polygon calibration mode, later inspect `presences_outside_calibration_area.csv` as well.
5. Record the starting and retained sample sizes. The full fitting stage requires at least 30 retained presence cells in this version; meeting this code minimum does not establish sampling adequacy.

**Check:** exclusions have been reviewed and any source corrections are documented. Retain original row identifiers and the unchanged original files.

### 4.6 Fit and evaluate SDMs locally

1. After accepting the preparation audit, change `RUN_MODE` to `"all"` and rerun in a fresh session with the same unchanged configuration and inputs. Preparation-to-fitting switch is supported within the same run directory.
2. Inspect background points, `outer_fold_counts.csv` and fold assignments. Verify adequate representation and inspect their geography; assess the provisional block size rather than treating it as measured spatial dependence.
3. Confirm predictor screening, encoding and hyperparameter selection use training partitions. The code implements correlation screening at 0.70 and a VIF limit of 5.
4. Review model and ensemble results using `spatial_CV_metrics_all_folds.csv`, `spatial_CV_summary.csv`, saved test predictions and held-out permutation importance.
5. Examine AUC, TSS, sensitivity, specificity, balanced accuracy, Boyce and prediction coverage together. Retain failed or missing metrics and diagnose their cause; do not select only favourable folds.
6. Examine final selected predictors, ensemble weights and the saved binary threshold under each species' `final_fit/` directory. Training-derived thresholds and weights must remain separate from outer test evaluation.

**Check:** every intended species completed successfully. Inspect `failed_species.csv` and species `FAILED.txt` files if present. A run that continues after a species error is not evidence of complete results.

### 4.7 Examine current and optional future maps

1. Inspect current continuous/binary suitability, algorithm disagreement, environmental novelty and area summaries.
2. Interpret continuous predictions as relative suitability. The novelty flag detects univariate range departures; it is not a calibrated confidence interval or a full multivariate novelty assessment.
3. For future scenarios, verify matching predictor identities, units, geometry and period/SSP/GCM metadata. Document every static non-climatic assumption; substitution is disabled by default.
4. Review per-GCM predictions and grouped scenario summaries separately. Do not mix SSPs or periods in a consensus.
5. For range maps, the code uses `current_binary + 2 * future_binary`: 0 = unsuitable in both, 1 = loss, 2 = gain, 3 = stable suitable. Check the legend against this encoding.
6. Retain suitable-area, GCM spread, suitable-vote fractions and agreement products with the scenario membership list.

**Check:** scenario metadata and map legends match the files. Geographic co-suitability is not proof of observed co-occurrence or malaria transmission.

### 4.8 Calculate and interpret overlap

**Nine-species outputs:** the SDM script produces area-weighted geographic Jaccard/Sørensen comparisons and optional two-axis PCA-env overlap. Preserve those methods as separate analyses; the PCA-env outputs are not the SVM metrics below.

**Five-species procedure:**

1. Identify the exact source SDM run for each species. Supply *An. gambiae* outputs separately and confirm all five runs use comparable domains, grids, periods and threshold provenance.
2. In a fresh session configure `environment_root`, `niche_out_dir`, `inputs` and `final_predictors` in the separate overlap script. Resolve predictor filename differences described in the data guide.
3. Use existing species-specific binary thresholds. Do not silently replace them with a generic 0.5 threshold. Supply single-band maps with identical geometry and valid 0/1 binary values.
4. Examine the common valid domain and occurrence audit. A changed domain changes areas and overlap estimates.
5. Run the script and inspect all ten rows of `pairwise_niche_overlap.csv`, including `SVM_status`, missing metrics and sample sizes.
6. Interpret the metrics using the table below. Retain `analysis_settings.csv`, standardization parameters, environmental centroids and `R_session_info.txt`.
7. Set `save_overlap_maps <- TRUE` only when map exports are required. The encoding is 0 = neither species suitable, 1 = species A only, 2 = species B only, 3 = both suitable; cells outside the common domain are NA.

| Metric | Implemented definition | Interpretation |
| --- | --- | --- |
| Geographic Sørensen | Twice shared suitable area / sum of both suitable areas | Higher values indicate more geographic overlap on the chosen common domain |
| SVM CV accuracy | Balanced species samples; radial SVM; same-cell grouped random folds | Lower accuracy indicates less separability for this classifier and predictor set; these are not spatial-block validation folds |
| Centroid distance | RMS difference between suitability-weighted standardized continuous-environment centroids | Lower values indicate closer environmental centres; the units are not geographic kilometres |

The retained script has SVM cutoff 0.80, unset Sørensen/centroid cutoffs and a two-of-three rule. These are exploratory choices. Continuous metrics can be reported without choosing new cutoffs. If categorical labels are retained, justify numerical rules independently, record them before the classified analysis and report sensitivity to those choices. Missing cutoffs yield unresolved labels; an unresolved label is not a non-overlap conclusion. There is no formal equivalency test in this SVM procedure.

The two scripts also differ in aspect handling: the overlap script's fallback aligns aspect angles before conversion. Establish compatibility with the original fitted SDM preprocessing, or supply documented north/east component rasters matching that run, before treating the workflows as one reproducible chain.

**Check:** each output is associated with the correct method, species list and source maps. Do not describe the nine-species and five-species overlap results as interchangeable.

### 4.9 Archive, review and release

1. Complete `metadata/run_log.csv` and `metadata/qa_checklist.md`. Record deviations, failed stages and decisions.
2. Preserve exact configured scripts, input checksums and locations, selected settings, package environment, source data references and complete outputs outside Git when necessary.
3. Link each manuscript figure to its run ID, source output and plotting/export method. Record any QGIS styling or manual layout operations separately.
4. Obtain technical review of the scientific settings and unresolved implementation differences. Attach reviewer and date to the record.
5. Commit code, documentation and appropriate metadata. Verify staged files before pushing. Keep generated rasters, model binaries and original occurrence data out of routine commits unless their distribution is intentionally arranged.
6. Confirm contributor attribution, dataset terms, code license and the intended manuscript citation. No license, DOI or institutional approval is supplied by this draft.

**Acceptance:** a reproduction is verified only after another clean local execution or reviewer reproduction successfully completes the intended stages and reproduces agreed outputs. This repository preparation has not performed that validation.

## 5. Troubleshooting and deviations

| Issue | Action |
| --- | --- |
| Missing raster or ambiguous filename | Compare stems with the correct script's input table; do not substitute a different product silently |
| Checkpoint fingerprint mismatch | Preserve the old run; choose a new output directory and document the change |
| Too few occurrences or inadequate folds | Inspect exclusions, sampling and validation design; record any justified redesign |
| Incompatible overlap maps | Return to input preparation and establish a consistent geometry/domain; do not suppress the check |
| Missing future predictors | Supply actual future inputs or explicitly document/enable a static-landscape assumption |
| Unclassified overlap | Report continuous metrics and missing criteria; do not invent thresholds to force a label |
| Package installation fails | Record R/platform details and resolve the dependency; do not replace the algorithm without documenting the change |

## 6. Documentation references

These references support project setup and software use. They do not establish universal biological cutoffs or verify the supplied datasets.

- [GitHub: creating a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository).
- [renv: environments, snapshot and restore](https://rstudio.github.io/renv/articles/renv.html).
- [CHELSA: dataset catalogue](https://www.chelsa-climate.org/). Select the citation and version belonging to the actual product used.
- [terra reference](https://rspatial.github.io/terra/reference/index.html), [maxnet manual](https://cran.r-project.org/web/packages/maxnet/maxnet.pdf), and [XGBoost R documentation](https://xgboost.readthedocs.io/en/stable/r_docs/).

## 7. Revision history

| Version | Date | Change |
| --- | --- | --- |
| 0.1 | 2026-10-29 | First SOP based on the two retrieved scripts; no model execution or technical approval |
