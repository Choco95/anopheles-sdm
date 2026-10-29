# Checks performed during repository preparation

Date: 29 October 2026.

- Both included analysis scripts were compared byte for byte with the retrieved originals; they match.
- Relative Markdown links were checked against included files.
- Metadata CSV row widths were checked against their headers.
- Git ignore rules were checked for occurrence inputs, raster inputs, model results, local logs and environment credentials.
- Staged files were checked for whitespace errors.
- The deliverable ZIP was checked for archive integrity and excludes Git history, original datasets and generated model results.

R is unavailable in the preparation environment. R parsing, package installation, input raster checks, model fitting, future projection and numerical reproduction were **not executed**. The supplied `scripts/00_check_environment.R` must be run locally, followed by the data audit and intended analyses. These preparation checks do not validate the scientific models.
