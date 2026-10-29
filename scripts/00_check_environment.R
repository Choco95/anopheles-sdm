# Parse without running: sourcing either analysis script would execute it.
scripts <- c("scripts/Anopheles_Africa_SDM_updated.R",
             "scripts/Anopheles_niche_overlap_5_species.R")
if (!all(file.exists(scripts))) stop("Open the project and run from its root folder.")
for (path in scripts) {
  invisible(parse(file = path))
  message("Parsed: ", path)
}
packages <- c("terra", "sf", "mgcv", "ranger", "maxnet", "xgboost", "pROC",
              "ecospat", "rnaturalearth", "rnaturalearthdata", "digest", "e1071")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(Package = packages, Available = available), row.names = FALSE)
if (!all(available)) stop("Install missing dependencies using scripts/00_setup.R.")
message("Code parsing and package availability checks passed. Input data and model execution remain unchecked.")
