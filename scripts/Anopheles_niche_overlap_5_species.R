# ANOPHELES NICHE OVERLAP: 5 SPECIES, 10 PAIRS
# Uses existing SDM maps; it does NOT refit the SDMs.
# Packages needed once: install.packages(c("terra", "e1071"))
# Run in RStudio. Environmental predictors are rebuilt from your raster folder.
# Results: CSV tables + 3 heatmaps. Optional geographic overlap GeoTIFFs.

# 1. SETTINGS ---------------------------------------------------------------
niche_out_dir <- file.path("E:/Anopheles/Anopheles_niche_overlap")
seed <- 42L
svm_folds <- 5L
svm_repeats <- 1L                 # increase to 3 for a sensitivity check
svm_max_per_species <- 500L       # balanced sample for each species pair
reference_sample_size <- 10000L   # one draw, then discard incomplete cells
svm_accuracy_cutoff <- 0.80
sorensen_cutoff <- NA_real_       # enter a justified, predefined cutoff
centroid_distance_cutoff <- NA_real_ # distance in standardized environmental units
save_overlap_maps <- FALSE       # TRUE saves one overlap map per pair
csv_crs <- "EPSG:4326"            # original WGS84 longitude/latitude coordinates
categorical_predictors <- "LULC"

# Read the selected predictors from your original environmental files.
# Keep your supplied species/map/occurrence paths below; no file dialogs.
environment_root <- "E:/Anopheles/Environmental data"
bioclim_folder <- file.path(environment_root, "bio_clim_var")
final_predictors <- c("bio2", "bio6", "bio14", "bio15", "Humidity",
  "Elevation", "Slope", "AspectEast", "AspectNorth", "TWI",
  "DistanceToWater", "LULC", "PopulationDensity", "HumanFootprint")

inputs <- data.frame(
  Species = c("arabiensis", "coluzzii", "coustani", "funestus", "gambiae"),
  Occurrences = c(
    "E:/Anopheles/Outputs_Arabiensis/Arabiensis/Arabiensis.csv",
    "E:/Anopheles/Outputs_Coluzzi/Coluzzi/Coluzzi.csv",
    "E:/Anopheles/Outputs_Coustani/Coustani/Coustani.csv",
    "E:/Anopheles/Outputs_funestus/Funestus/Funestus.csv",
    "E:/Anopheles/Outputs_gambiae/gambiae/Gambiae.csv"
  ),
  Continuous =  c(
    "E:/Anopheles/Outputs_Arabiensis/Arabiensis/ensemble_map.tif",
    "E:/Anopheles/Outputs_Coluzzi/Coluzzi/ensemble_map.tif",
    "E:/Anopheles/Outputs_Coustani/Coustani/ensemble_map.tif",
    "E:/Anopheles/Outputs_funestus/Funestus/ensemble_map.tif",
    "E:/Anopheles/Outputs_gambiae/gambiae/ensemble_map.tif"
  ), # ensemble_map.tif
  Binary =  c(
    "E:/Anopheles/Outputs_Arabiensis/Arabiensis/binary_map.tif",
    "E:/Anopheles/Outputs_Coluzzi/Coluzzi/binary_map.tif",
    "E:/Anopheles/Outputs_Coustani/Coustani/binary_map.tif",
    "E:/Anopheles/Outputs_funestus/Funestus/binary_map.tif",
    "E:/Anopheles/Outputs_gambiae/gambiae/binary_map.tif"
  ),     # binary_map.tif (0/1), existing SDM thresholds
  stringsAsFactors = FALSE
)
# Use full original study-area maps for every species, with matching masks.
# Do not mix filtered and unfiltered maps or use a new threshold of 0.5.

# 2. PACKAGES AND INPUT FILES ------------------------------------------------
needed <- c("terra", "e1071")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages first: install.packages(c(",
                         paste(sprintf('"%s"', missing), collapse = ", "), "))")
dir.create(niche_out_dir, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(progress = 1)
if (svm_folds < 2 || svm_repeats < 1 || svm_max_per_species < 10) stop("Invalid SVM fold/repeat/sample settings.")
if (!is.finite(svm_accuracy_cutoff) || svm_accuracy_cutoff < 0 || svm_accuracy_cutoff > 1) stop("SVM cutoff must be in [0,1].")
if (is.finite(sorensen_cutoff) && (sorensen_cutoff < 0 || sorensen_cutoff > 1)) stop("Sorensen cutoff must be in [0,1].")
if (is.finite(centroid_distance_cutoff) && centroid_distance_cutoff < 0) stop("Centroid cutoff must be nonnegative.")
start_time <- Sys.time()
set.seed(seed)
if (any(!file.exists(unlist(inputs[c("Occurrences", "Continuous", "Binary")])))) {
  stop("One or more occurrence/map paths do not exist.")
}
# BUILD THE SHARED ENVIRONMENTAL STACK FROM THE ORIGINAL FILES
# The first species' ensemble supplies the exact grid; no Africa shapefile
# or existing env_sdm object is needed. The later common-domain step masks
# the comparison to valid cells shared by all species and all predictors.
if (!dir.exists(environment_root)) stop("Environmental folder not found: ", environment_root)
if (!dir.exists(bioclim_folder)) stop("Bioclim folder not found: ", bioclim_folder)
reference_grid <- terra::rast(inputs$Continuous[1])
if (terra::nlyr(reference_grid) != 1L || !nzchar(terra::crs(reference_grid))) {
  stop("The first ensemble map must be a single-band raster with a defined CRS.")
}
aligned_dir <- file.path(niche_out_dir, "environment_aligned")
dir.create(aligned_dir, recursive = TRUE, showWarnings = FALSE)
root_files <- list.files(environment_root, pattern = "\\.tif(f)?$",
                         full.names = TRUE, ignore.case = TRUE)
find_environment_file <- function(variable, required = TRUE) {
  stems <- tools::file_path_sans_ext(basename(root_files))
  hits <- root_files[tolower(stems) == tolower(variable)]
  if (length(hits) > 1L) stop("More than one raster found for ", variable, ": ", paste(hits, collapse = ", "))
  if (!length(hits)) {
    if (required) stop("Missing raster: ", file.path(environment_root, paste0(variable, ".tif")))
    return(NULL)
  }
  hits[1]
}
align_environment <- function(r, variable, method = "bilinear") {
  if (terra::nlyr(r) != 1L) stop(variable, " must be a single-band raster.")
  if (!nzchar(terra::crs(r))) stop("Original CRS missing from predictor: ", variable)
  if (!terra::compareGeom(reference_grid, r, stopOnError = FALSE)) {
    message("Aligning ", variable, " to the ensemble grid...")
    r <- terra::project(r, reference_grid, method = method,
      filename = file.path(aligned_dir, paste0(variable, ".tif")), overwrite = TRUE)
  }
  names(r) <- variable
  r
}
message("Loading the four retained bioclimatic predictors...")
bio_files <- list.files(bioclim_folder, pattern = "\\.tif(f)?$",
                        full.names = TRUE, ignore.case = TRUE)
if (!length(bio_files)) stop("No bioclimatic GeoTIFFs found in ", bioclim_folder)
bio_labels <- tolower(tools::file_path_sans_ext(basename(bio_files)))
bio_numbers <- suppressWarnings(as.integer(sub(".*bio[_ -]?0*([0-9]+).*", "\\1", bio_labels)))
required_bio <- c(2L, 6L, 14L, 15L)
# Most original SDM folders have one file per bioclimatic variable.
# Also accept one multiband file if its layer names explicitly identify bioN.
if (length(bio_files) == 1L && terra::nlyr(terra::rast(bio_files)) > 1L) {
  bio_all <- terra::rast(bio_files)
  layer_labels <- tolower(names(bio_all))
  bio_numbers <- suppressWarnings(as.integer(sub(".*bio[_ -]?0*([0-9]+).*", "\\1", layer_labels)))
  selected_bio <- lapply(required_bio, function(n) {
    k <- which(bio_numbers == n)
    if (length(k) != 1L) stop("Cannot identify bio", n, " uniquely in the multiband bioclim raster. Set its layer names explicitly.")
    align_environment(bio_all[[k]], paste0("bio", n))
  })
} else {
  selected_bio <- lapply(required_bio, function(n) {
    k <- which(bio_numbers == n)
    if (length(k) != 1L) stop("Cannot identify one file for bio", n, ". Files: ", paste(basename(bio_files), collapse = ", "))
    align_environment(terra::rast(bio_files[k]), paste0("bio", n))
  })
}
names(selected_bio) <- paste0("bio", required_bio)
non_bio <- c("Humidity", "Elevation", "Slope", "TWI", "DistanceToWater",
             "LULC", "PopulationDensity", "HumanFootprint")
selected_other <- lapply(non_bio, function(v) {
  align_environment(terra::rast(find_environment_file(v)), v,
                    method = if (v == "LULC") "near" else "bilinear")
})
names(selected_other) <- non_bio
# Use existing aspect-component rasters if both are available. Otherwise
# reproduce your original SDM conversion: aspect in DEGREES -> cos/sin,
# after aligning Aspect.tif to the same grid as in the SDM workflow.
north_path <- find_environment_file("AspectNorth", required = FALSE)
east_path <- find_environment_file("AspectEast", required = FALSE)
if (!is.null(north_path) && !is.null(east_path)) {
  aspect_north <- align_environment(terra::rast(north_path), "AspectNorth")
  aspect_east <- align_environment(terra::rast(east_path), "AspectEast")
} else {
  aspect <- align_environment(terra::rast(find_environment_file("Aspect")), "Aspect")
  aspect_north <- cos(aspect * pi / 180)
  aspect_east <- sin(aspect * pi / 180)
  names(aspect_north) <- "AspectNorth"; names(aspect_east) <- "AspectEast"
}
# c() combines the rasters in a known order; each named layer is selected
# explicitly to avoid accidentally including wind, NDVI or other variables.
env <- c(terra::rast(selected_bio), terra::rast(selected_other),
         aspect_north, aspect_east)
env <- env[[final_predictors]]
if (anyDuplicated(names(env))) stop("Predictor names must be unique.")
message("Shared environmental stack built: ", terra::nlyr(env), " predictors.")
print(names(env))

write.csv(inputs, file.path(niche_out_dir, "input_files.csv"), row.names = FALSE)
write.csv(data.frame(Predictor = names(env),
                     Categorical = names(env) %in% categorical_predictors),
          file.path(niche_out_dir, "predictors_used.csv"), row.names = FALSE)

# 3. CHECK MAPS AND DEFINE ONE COMMON VALID DOMAIN ---------------------------
# Geometry is checked rather than silently resampling masks or binary maps.
# Multi-band map inputs are rejected to avoid confusing suitability with SD.
read_map <- function(path, kind) {
  r <- terra::rast(path)
  if (terra::nlyr(r) != 1L) stop(kind, " must be a single-band raster: ", path)
  if (!nzchar(terra::crs(r))) stop("Missing map CRS: ", path)
  if (!terra::compareGeom(env[[1]], r, stopOnError = FALSE)) {
    stop("Map grid/CRS/extent differs from the rebuilt environmental stack: ", path,
         ". Align the inputs on the same grid first.")
  }
  limits <- terra::global(r, c("min", "max"), na.rm = TRUE)
  if (any(!is.finite(as.numeric(unlist(limits[1, ], use.names = FALSE))))) stop("Map has no valid values: ", path)
  if (limits$min < -1e-7 || limits$max > 1 + 1e-7) {
    stop(kind, " must have values in [0, 1]: ", path)
  }
  if (kind == "Binary") {
    bad <- terra::global((r != 0) & (r != 1), "sum", na.rm = TRUE)[1, 1]
    if (bad > 0) stop("Binary raster contains values other than 0/1: ", path)
  }
  r
}
message("Checking the ten map inputs...")
continuous <- lapply(inputs$Continuous, read_map, kind = "Continuous")
binary <- lapply(inputs$Binary, read_map, kind = "Binary")
names(continuous) <- names(binary) <- inputs$Species
cont_stack <- terra::rast(continuous); names(cont_stack) <- paste0("pred_", inputs$Species)
bin_stack <- terra::rast(binary); names(bin_stack) <- paste0("bin_", inputs$Species)
message("Preparing the common valid domain; this reads the rasters once.")
all_inputs <- c(env, cont_stack, bin_stack)
missing_count <- terra::app(is.na(all_inputs), "sum")
common_mask <- terra::ifel(missing_count == 0, 1, NA,
  filename = file.path(niche_out_dir, "common_valid_domain.tif"), overwrite = TRUE)
area <- terra::mask(terra::cellSize(env[[1]], unit = "km"), common_mask)
common_area_km2 <- terra::global(area, "sum", na.rm = TRUE)[1, 1]
if (!is.finite(common_area_km2) || common_area_km2 <= 0) stop("No common valid domain.")
message("Common area: ", round(common_area_km2), " km2")

# 4. EXTRACT PRESENCE ENVIRONMENTS -------------------------------------------
# Fresh CSVs avoid the earlier SpatialPointsDataFrame/geom mismatch.
# No record is modified in the source files. Keep one record per environmental
# cell for this comparison; report all exclusions in occurrence_audit.csv.
read_presence <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(gsub("[^a-z0-9]", "", tolower(names(d))))
  find_col <- function(candidates) {
    match_id <- match(candidates, nm); match_id <- match_id[!is.na(match_id)]
    if (!length(match_id)) stop("Coordinate columns not detected: ", path)
    match_id[1]
  }
  x <- suppressWarnings(as.numeric(as.character(d[[find_col(c("longitude", "decimallongitude", "lon", "long", "x"))]])))
  y <- suppressWarnings(as.numeric(as.character(d[[find_col(c("latitude", "decimallatitude", "lat", "y"))]])))
  valid <- is.finite(x) & is.finite(y)
  if (csv_crs == "EPSG:4326") valid <- valid & abs(x) <= 180 & abs(y) <= 90
  if (!any(valid)) stop("No valid coordinates: ", path)
  p <- terra::vect(data.frame(x = x[valid], y = y[valid]), geom = c("x", "y"), crs = csv_crs)
  p <- terra::project(p, terra::crs(env))
  cell <- terra::cellFromXY(env[[1]], terra::crds(p))
  e <- as.data.frame(terra::extract(env, p, ID = FALSE, raw = TRUE))
  mask_value <- terra::extract(common_mask, p, ID = FALSE)[[1]]
  ok <- !is.na(cell) & is.finite(mask_value) & complete.cases(e)
  ok <- ok & apply(as.matrix(e), 1, function(z) all(is.finite(z)))
  retained <- which(ok); retained <- retained[!duplicated(cell[retained])]
  if (length(retained) < 10L) stop("Fewer than 10 valid distinct-cell presences: ", path)
  list(env = e[retained, , drop = FALSE], cell = cell[retained],
       audit = data.frame(N_input = nrow(d), N_valid_coordinates = sum(valid),
         N_on_common_domain = sum(ok), N_unique_cells = length(retained)))
}
presences <- lapply(inputs$Occurrences, read_presence)
names(presences) <- inputs$Species
occ_audit <- do.call(rbind, lapply(seq_along(presences), function(i) {
  cbind(Species = inputs$Species[i], presences[[i]]$audit)
}))
write.csv(occ_audit, file.path(niche_out_dir, "occurrence_audit.csv"), row.names = FALSE)

# 5. SUITABILITY-WEIGHTED ENVIRONMENTAL CENTROIDS -----------------------------
# A SINGLE common random draw for all species; no repeated search for NA cells.
# terra's geographic random sampling accounts for latitude-based cell area.
# Centroids are in standardized environmental space, NOT geographic km.
# LULC is categorical: included as dummy variables in the SVM but excluded
# from the numerical environmental centroid distance.
message("Sampling the common environmental space for centroids...")
set.seed(seed)
# Cell-only sampling can return a vector or a one-column matrix.
# Convert directly to numeric instead of assuming a named data-frame column.
sampled_cells <- as.numeric(terra::spatSample(env[[1]],
  size = reference_sample_size, method = "random", na.rm = FALSE,
  cells = TRUE, values = FALSE, as.df = FALSE))
# raw=TRUE keeps land-cover category codes numeric during extraction.
reference <- as.data.frame(terra::extract(c(env, cont_stack, common_mask),
  sampled_cells, raw = TRUE))
reference <- reference[complete.cases(reference), , drop = FALSE]
reference <- reference[apply(as.matrix(reference), 1, function(z) all(is.finite(z))), , drop = FALSE]
if (nrow(reference) < 100L) stop("Fewer than 100 valid reference samples. Crop to the common study area or increase reference_sample_size.")
numeric_predictors <- setdiff(names(env), categorical_predictors)
ref_env <- as.matrix(reference[, numeric_predictors, drop = FALSE])
ref_sd <- apply(ref_env, 2, sd)
keep <- is.finite(ref_sd) & ref_sd > 0
if (!any(keep)) stop("No varying continuous predictors for centroid distances.")
ref_env <- ref_env[, keep, drop = FALSE]
Z <- scale(ref_env)
centroids <- t(vapply(inputs$Species, function(sp) {
  w <- reference[[paste0("pred_", sp)]]
  if (sum(w) <= 0) stop("Zero suitability weights for ", sp)
  colSums(Z * w) / sum(w)
}, numeric(ncol(Z))))
rownames(centroids) <- inputs$Species; colnames(centroids) <- colnames(Z)
write.csv(data.frame(Species = rownames(centroids), centroids, check.names = FALSE),
          file.path(niche_out_dir, "environmental_centroids.csv"), row.names = FALSE)
write.csv(data.frame(Predictor = colnames(Z), Centre = attr(Z, "scaled:center"),
                     Scale_SD = attr(Z, "scaled:scale")),
          file.path(niche_out_dir, "centroid_standardization.csv"), row.names = FALSE)

# 6. CROSS-VALIDATED SVM SEPARABILITY ----------------------------------------
# Equal species sample sizes; fixed radial kernel and cost = 1, no tuning.
# Predictor scaling is learned INSIDE each training fold.
# The same raster cell is always assigned to the same fold across both species.
# These are grouped random folds, not spatial-block validation.
# Low accuracy indicates limited separability for this classifier/settings;
# it is evidence to combine with other metrics, not a formal overlap test.
svm_pair <- function(a, b, pair_seed) {
  set.seed(pair_seed)
  n <- min(svm_max_per_species, nrow(a$env), nrow(b$env))
  ia <- sample(seq_len(nrow(a$env)), n); ib <- sample(seq_len(nrow(b$env)), n)
  d <- rbind(a$env[ia, , drop = FALSE], b$env[ib, , drop = FALSE])
  cells <- c(a$cell[ia], b$cell[ib]); y <- factor(rep(c("A", "B"), each = n))
  for (v in categorical_predictors) d[[v]] <- factor(d[[v]])
  cats <- categorical_predictors[vapply(d[categorical_predictors], nlevels, integer(1)) < 2L]
  if (length(cats)) d <- d[, setdiff(names(d), cats), drop = FALSE]
  X <- model.matrix(~ . - 1, data = d)
  if (!ncol(X)) stop("No predictors for SVM.")
  unique_cells <- unique(cells)
  kind <- vapply(unique_cells, function(cl) paste(sort(unique(as.character(y[cells == cl]))), collapse = ""), character(1))
  k <- min(svm_folds, n)
  accuracy <- balanced <- numeric(svm_repeats)
  for (rep_id in seq_len(svm_repeats)) {
    cell_fold <- integer(length(unique_cells))
    for (stratum in unique(kind)) {
      ids <- which(kind == stratum)
      cell_fold[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
    }
    fold <- cell_fold[match(cells, unique_cells)]
    oof <- rep(NA_character_, length(y))
    for (f in sort(unique(fold))) {
      train <- which(fold != f); test <- which(fold == f)
      if (length(unique(y[train])) < 2L) stop("A training fold has only one species.")
      mu <- colMeans(X[train, , drop = FALSE]); ss <- apply(X[train, , drop = FALSE], 2, sd)
      active <- is.finite(ss) & ss > 0
      if (!any(active)) stop("No varying predictors in an SVM training fold.")
      xtrain <- sweep(sweep(X[train, active, drop = FALSE], 2, mu[active], "-"), 2, ss[active], "/")
      xtest <- sweep(sweep(X[test, active, drop = FALSE], 2, mu[active], "-"), 2, ss[active], "/")
      m <- e1071::svm(x = xtrain, y = y[train], type = "C-classification",
        kernel = "radial", cost = 1, gamma = 1 / ncol(xtrain), scale = FALSE,
        probability = FALSE)
      oof[test] <- as.character(predict(m, xtest))
    }
    if (anyNA(oof)) stop("Incomplete out-of-fold SVM predictions.")
    accuracy[rep_id] <- mean(oof == as.character(y))
    balanced[rep_id] <- mean(vapply(levels(y), function(cl) mean(oof[y == cl] == cl), numeric(1)))
  }
  list(accuracy = mean(accuracy), balanced_accuracy = mean(balanced), n = n,
       sd = if (svm_repeats > 1L) sd(accuracy) else NA_real_)
}

# 7. TEN PAIRS: SORENSEN, SVM, CENTROID DISTANCE -----------------------------
# Geographic Sørensen uses cell AREA, not raw pixel counts, on the same domain:
# 2 * intersecting suitable area / (species-A suitable area + species-B area).
# Environmental distance = root-mean-square difference across standardized
# centroid coordinates. Lower values indicate more similar niche positions.
message("Calculating suitable areas...")
suitable_area <- vapply(binary, function(r) terra::global(area * r, "sum", na.rm = TRUE)[1, 1], numeric(1))
pairs <- combn(seq_len(nrow(inputs)), 2)
results <- vector("list", ncol(pairs))
for (p in seq_len(ncol(pairs))) {
  i <- pairs[1, p]; j <- pairs[2, p]; a <- inputs$Species[i]; b <- inputs$Species[j]
  tick <- Sys.time()
  message("Pair ", p, "/", ncol(pairs), ": ", a, " vs ", b)
  overlap_area <- terra::global(area * binary[[i]] * binary[[j]], "sum", na.rm = TRUE)[1, 1]
  denominator <- suitable_area[i] + suitable_area[j]
  sorensen <- if (denominator > 0) 2 * overlap_area / denominator else NA_real_
  centroid_distance <- sqrt(mean((centroids[i, ] - centroids[j, ])^2))
  svm_result <- tryCatch(svm_pair(presences[[i]], presences[[j]], seed + p),
    error = function(e) list(accuracy = NA_real_, balanced_accuracy = NA_real_,
      n = NA_integer_, sd = NA_real_, error = conditionMessage(e)))
  c_svm <- if (is.finite(svm_result$accuracy)) svm_result$accuracy <= svm_accuracy_cutoff else NA
  c_sor <- if (is.finite(sorensen) && is.finite(sorensen_cutoff)) sorensen >= sorensen_cutoff else NA
  c_cent <- if (is.finite(centroid_distance) && is.finite(centroid_distance_cutoff)) centroid_distance <= centroid_distance_cutoff else NA
  criteria <- c(c_svm, c_sor, c_cent)
  n_true <- sum(criteria, na.rm = TRUE); n_unknown <- sum(is.na(criteria))
  classification <- if (n_true >= 2L) "Overlap supported" else if (n_true + n_unknown < 2L) {
    "Insufficient support"
  } else "Not classified: missing cutoff or metric"
  results[[p]] <- data.frame(Species_A = a, Species_B = b,
    SVM_CV_accuracy = svm_result$accuracy, SVM_balanced_accuracy = svm_result$balanced_accuracy,
    SVM_repeat_SD = svm_result$sd, SVM_n_per_species = svm_result$n,
    Sorensen_area = sorensen, Overlap_km2 = overlap_area,
    Suitable_A_km2 = suitable_area[i], Suitable_B_km2 = suitable_area[j],
    Centroid_distance_RMS = centroid_distance, Criterion_SVM = c_svm,
    Criterion_Sorensen = c_sor, Criterion_centroid = c_cent,
    Criteria_met = n_true, Criteria_unknown = n_unknown,
    Overlap_classification = classification,
    SVM_status = if (is.null(svm_result$error)) "OK" else svm_result$error,
    Seconds = as.numeric(difftime(Sys.time(), tick, units = "secs")),
    stringsAsFactors = FALSE, row.names = NULL)
  # Save after every pair so completed work is retained if interrupted.
  write.csv(do.call(rbind, results[seq_len(p)]), file.path(niche_out_dir, "pairwise_niche_overlap.csv"), row.names = FALSE)
  if (save_overlap_maps) {
    # 0 = neither suitable, 1 = A only, 2 = B only, 3 = both suitable.
    r <- terra::mask(binary[[i]] + 2 * binary[[j]], common_mask)
    terra::writeRaster(r, file.path(niche_out_dir, paste0(a, "_", b, "_overlap.tif")), overwrite = TRUE)
  }
  message(round(100 * p / ncol(pairs)), "% of pairs completed")
}
results <- do.call(rbind, results)

# 8. SUMMARY FIGURES AND SETTINGS -------------------------------------------
heatmap_metric <- function(column, title, file, diagonal, bounds = NULL) {
  M <- matrix(NA_real_, nrow(inputs), nrow(inputs), dimnames = list(inputs$Species, inputs$Species))
  diag(M) <- diagonal
  for (p in seq_len(nrow(results))) {
    i <- match(results$Species_A[p], inputs$Species); j <- match(results$Species_B[p], inputs$Species)
    M[i, j] <- M[j, i] <- results[[column]][p]
  }
  write.csv(M, file.path(niche_out_dir, paste0(file, ".csv")), row.names = TRUE)
  png(file.path(niche_out_dir, paste0(file, ".png")), width = 2100, height = 1800, res = 250)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(8, 8, 4, 2))
  pal <- colorRampPalette(c("#F1F7F8", "#82C4BB", "#087B75", "#064B48"))(100)
  limits <- if (is.null(bounds)) range(M, na.rm = TRUE) else bounds
  if (diff(limits) == 0) limits <- limits + c(-0.01, 0.01)
  image(seq_len(nrow(M)), seq_len(ncol(M)), t(M[nrow(M):1, , drop = FALSE]),
    col = pal, zlim = limits, axes = FALSE, xlab = "", ylab = "", main = title)
  axis(1, at = seq_len(nrow(M)), labels = inputs$Species, las = 2, tick = FALSE)
  axis(2, at = seq_len(nrow(M)), labels = rev(inputs$Species), las = 2, tick = FALSE)
  for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M))) {
    value <- M[i, j]
    text(j, nrow(M) + 1 - i, if (is.finite(value)) sprintf("%.3f", value) else "—",
      col = if (is.finite(value) && value > mean(limits)) "white" else "#19384A", cex = 1.1)
  }
  mtext("Anopheles species; numbers are the pairwise metric values", side = 1, line = 6, cex = .75)
}
heatmap_metric("SVM_CV_accuracy", "SVM separability: lower accuracy = less separation", "SVM_accuracy_heatmap", NA_real_, c(0, 1))
heatmap_metric("Sorensen_area", "Geographic Sørensen similarity: higher = more overlap", "Sorensen_heatmap", 1, c(0, 1))
heatmap_metric("Centroid_distance_RMS", "Environmental centroid distance: lower = closer", "Centroid_distance_heatmap", 0)
settings <- data.frame(
  Setting = c("Seed", "SVM_folds", "SVM_repeats", "SVM_max_per_species", "Reference_sample_requested", "Reference_sample_valid", "Common_area_km2", "SVM_accuracy_cutoff", "Sorensen_cutoff", "Centroid_distance_cutoff", "Centroid_dimensions", "Centroid_definition", "SVM_definition"),
  Value = as.character(c(seed, svm_folds, svm_repeats, svm_max_per_species,
    reference_sample_size, nrow(reference), common_area_km2, svm_accuracy_cutoff,
    sorensen_cutoff, centroid_distance_cutoff, ncol(Z),
    "Suitability-weighted, standardized continuous environments; RMS distance",
    "Balanced presences; same-cell grouped random CV; radial kernel; cost=1")))
write.csv(settings, file.path(niche_out_dir, "analysis_settings.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(niche_out_dir, "R_session_info.txt"))
print(results[, c("Species_A", "Species_B", "SVM_CV_accuracy", "Sorensen_area", "Centroid_distance_RMS", "Overlap_classification")])
message("Done in ", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), " minutes. Outputs: ", niche_out_dir)
if (!is.finite(sorensen_cutoff) || !is.finite(centroid_distance_cutoff)) {
  message("Metrics calculated. The 2-of-3 rule remains unresolved where numerical cutoffs are missing; no cutoffs have been invented.")
}
# Sources: https://search.r-project.org/CRAN/refmans/e1071/html/svm.html
# https://rspatial.github.io/terra/reference/sample.html
# https://rspatial.github.io/terra/reference/cellSize.html
# This is a transparent exploratory three-metric comparison, not a calibrated
# significance test or proof of niche equivalence. Random CV can retain spatial
# dependence; occurrence survey bias and SDM calibration affect the results.
