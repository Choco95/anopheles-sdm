# Run from the project root. Installs only packages that are currently missing.
# Once a validated renv.lock exists, prefer renv::restore() for replication.
packages <- c("terra", "sf", "mgcv", "ranger", "maxnet", "xgboost", "pROC",
              "ecospat", "rnaturalearth", "rnaturalearthdata", "digest", "e1071")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Packages still unavailable: ", paste(still_missing, collapse = ", "))
message("All declared analysis packages are available. Record versions after a successful local run.")
