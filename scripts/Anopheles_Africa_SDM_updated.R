# =============================================================================
# ANOPHELES AFRICA SDM — nine species, CHELSA working grid ~5 km
# Version 1.0 | 2026-10-02
# Prepared for LOCAL execution by the user. No models were run during preparation.
#
# HOW TO RUN
# 1. Edit ROOT, ENV_DIR, OCC_DIR and OUT_DIR below.
# 2. Put the nine original CSVs in OCC_DIR. Keep their original filenames.
# 3. Put CHELSA single-band BIO rasters in ENV_DIR/bio_clim_var and the other
#    rasters directly in ENV_DIR, with filenames matching the screenshot.
# 4. Run this script in RStudio. Missing packages produce an installation command.
# 5. First use RUN_MODE='prepare' to inspect the data/raster audit on YOUR machine.
#    Then set RUN_MODE='all' to fit, validate and map all nine species locally.
#    RUN_MODE is a workflow switch, not an invitation for an external agent to fit.
#
# INCLUDED
# - Explicit species-label filtering; original labels and excluded rows retained.
# - Africa mask; exact duplicate and one-per-grid-cell filtering within species.
# - Continuous/categorical alignment; aspect -> northness/eastness BEFORE warp.
# - Spatial outer CV; spatial inner tuning; fold-specific correlation/VIF selection.
# - Maxnet (MaxEnt without Java), ranger RF, XGBoost and mgcv GAM.
# - Training-only ensemble weights/thresholds; untouched outer-fold evaluation.
# - AUC, TSS, sensitivity, specificity, balanced accuracy, Boyce and coverage.
# - Held-out permutation importance; descriptive partial-dependence curves.
# - Disk-based raster prediction, model disagreement and extrapolation flags.
# - Optional per-GCM future projections, consensus, agreement and four-class shifts.
# - Area-weighted geographic overlap, richness, and common PCA-env niche overlap.
# - Species checkpoints; input/config/code fingerprints prevent stale restarts.
#
# SCIENTIFIC SETTINGS TO DOCUMENT
# - Grid: the existing CHELSA file geometry is used, not an invented resolution.
# - CALIBRATION_MODE='africa' is an explicit broad-domain baseline. For a justified
#   accessible area, choose 'polygons' and supply one polygon per species below.
# - Background defaults to random available cells, NOT confirmed absences.
#   Target-group sampling is optional and requires defensible common survey effort.
# - BLOCK_KM=100 is a starting design setting, NOT an empirically estimated range.
#   Inspect fold maps/counts; rerun with alternative sizes in separate output runs.
# - Class balancing changes score calibration. Outputs are relative suitability,
#   not mosquito abundance or calibrated malaria-transmission probabilities.
# - Thresholds optimized on background data depend on background definition.
# - Default taxonomy is strict. Complex/like labels do not identify sensu-stricto
#   species. Quadriannulatus sp. B is excluded (now An. amharicus).
# - Future non-climatic predictors may be held constant ONLY if explicitly enabled.
# - PCA-env is descriptive, occurrence-based, continuous-environment overlap in
#   two axes; it is not a niche-equivalency test or an SDM performance measure.
#
# REFERENCES / API DOCUMENTATION
# Valavi et al. 2019: https://doi.org/10.1111/2041-210X.13107
# Coetzee et al. 2013: https://doi.org/10.11646/zootaxa.3619.3.2
# Broennimann et al. 2012: https://doi.org/10.1111/j.1466-8238.2011.00698.x
# https://rspatial.github.io/terra/reference/extract.html
# https://rspatial.github.io/terra/reference/resample.html
# https://cran.r-project.org/web/packages/maxnet/maxnet.pdf
# https://xgboost.readthedocs.io/en/stable/r_docs/R-package/docs/reference/xgb.train.html
# =============================================================================

# 0. USER SETTINGS -------------------------------------------------------------
ROOT <- "D:/Anopheles"                    # CHANGE THIS
ENV_DIR <- file.path(ROOT, "Environmental_data")
OCC_DIR <- file.path(ROOT, "Occurrence_data")
OUT_DIR <- file.path(ROOT, "Outputs", "Africa_5km_v1")
RUN_MODE <- "prepare"                     # 'prepare' or 'all'

CFG <- list(
  seed = 20261002L,
  algorithms = c("MAXNET", "RF", "XGB", "GAM"),
  outer_k = 5L, inner_k = 3L,
  block_km = 100,                         # provisional; assess sensitivity
  min_test_presences = 5L, min_test_background = 20L,
  min_species_presences = 30L,
  background_n = 10000L,
  background_method = "random",          # 'random' or 'target_group'
  target_group_effort_confirmed = FALSE,
  calibration_mode = "africa",           # 'africa' or 'polygons'
  correlation_limit = 0.70, vif_limit = 5,
  selection_sample_n = 10000L,
  ensemble_weight = "TSS",               # 'TSS' or 'AUC_above_0.5'
  minimum_inner_coverage = 0.95,
  threads = 2L,                          # parallelize within RF/XGB only
  permutation_repeats = 3L,
  prediction_memfrac = 0.45,
  make_response_curves = TRUE,
  make_niche_overlap = TRUE,
  background_pca_n = 20000L,
  niche_grid_resolution = 100L,
  aspect_units = "degrees",              # confirm from raster metadata
  aspect_flat_code = -1,                 # only this sentinel is treated as flat
  water_nonwater_is_na = FALSE,           # TRUE only if NA explicitly means land
  # Optional local Africa boundary. Empty -> packaged Natural Earth 1:110m.
  # A more detailed local boundary is preferable for coastal/island records.
  africa_boundary = "",
  # No silent future substitution; set TRUE for an explicit static-landscape run.
  allow_static_future_nonbioclim = FALSE
)

SPECIES <- data.frame(
  species = c("funestus", "arabiensis", "coluzzii", "coustani", "merus",
              "quadriannulatus", "stephensi", "moucheti", "nili"),
  file = c("funestus.csv", "Arabiensis.csv", "Coluzzi.csv", "Coustani.csv",
           "Merus.csv", "quadra.csv", "Stephensi.csv", "moucheti.csv", "nili.csv"),
  stringsAsFactors = FALSE
)
# Exact accepted labels AFTER removal of an optional 'Anopheles ' prefix.
# Do not add complex/like labels unless the modelling unit is explicitly changed.
ACCEPTED_LABELS <- list(
  funestus = c("funestus", "funestus s.s."),
  arabiensis = "arabiensis",
  coluzzii = c("coluzzii", "coluzzi", "coluzzii_gambiae_m form"),
  coustani = c("coustani", "coustani s.s."),
  merus = "merus",
  quadriannulatus = c("quadriannulatus", "quadriannulatus sp. a"),
  stephensi = "stephensi", moucheti = "moucheti", nili = "nili"
)
CALIBRATION_FILES <- setNames(rep("", nrow(SPECIES)), SPECIES$species)
# For calibration_mode='polygons', supply each species' .gpkg/.shp path above.

# Predictor names are standardized here; the filenames retain the user's spelling.
EXTRA <- data.frame(
  name = c("lulc", "human_footprint", "humidity", "ndvi", "population_density",
           "slope", "solar_radiation", "elevation", "wind_speed", "aspect",
           "distance_to_water", "twi", "water"),
  stem = c("ESA_LULC_Africa", "Human_Footprint", "Humidity", "NDVI",
           "populationDensity", "Slope", "Solar Radiation", "wc2.1_5m_elev",
           "WindSpeed", "Aspect", "DistanceToWater", "TWI", "WaterRasterized"),
  type = c("categorical", rep("continuous", 8), "angle", "continuous",
           "continuous", "binary"),
  stringsAsFactors = FALSE
)
# Optional scaling applied AFTER reading (terra already honors supported GDAL
# scale/offset metadata). Do not apply a scale twice. Defaults do not guess units.
# Example ONLY if your file documentation requires it: INPUT_MULTIPLIER <- c(ndvi=0.0001)
INPUT_MULTIPLIER <- numeric(0)
INPUT_OFFSET <- numeric(0)

# Optional future rows; do not mix periods/SSPs in one consensus group.
# Each folder must have its own bio_clim_var subfolder containing BIO rasters.
# Supplied non-bioclim layers use the same stems as EXTRA. Missing ones are allowed
# only when CFG$allow_static_future_nonbioclim is TRUE; each substitution is logged.
FUTURE <- data.frame(gcm=character(), ssp=character(), period=character(),
                     folder=character(), stringsAsFactors=FALSE)
# Example (replace with real paths and your actual GCM names):
# FUTURE <- data.frame(
#   gcm=c('GCM_A','GCM_B','GCM_A','GCM_B'),
#   ssp=c('SSP245','SSP245','SSP585','SSP585'),
#   period=rep('2041_2070',4),
#   folder=file.path(ROOT,'Future',c('A_245','B_245','A_585','B_585')))

# 1. PACKAGES AND GENERAL HELPERS ----------------------------------------------
required <- c("terra", "sf", "mgcv", "ranger", "maxnet", "xgboost", "pROC",
              "ecospat", "rnaturalearth", "rnaturalearthdata", "digest")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop("Install missing packages first:\ninstall.packages(c(",
  paste(sprintf('"%s"', missing_pkgs), collapse=", "), "))", call.=FALSE)
if (!RUN_MODE %in% c("prepare", "all")) stop("Invalid RUN_MODE")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"temp"), showWarnings=FALSE)
terra::terraOptions(memfrac=CFG$prediction_memfrac,
                    tempdir=file.path(OUT_DIR,"temp"), progress=1)
AEA <- "+proj=aea +lat_1=-18 +lat_2=21 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m +no_defs"
WOPT <- list(datatype="FLT4S", gdal=c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
csv <- function(x, path) utils::write.csv(x, path, row.names=FALSE, na="")
safe_id <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)
finite_rows <- function(x) {
  if (!ncol(x)) return(rep(FALSE,nrow(x)))
  apply(as.matrix(x),1,function(z) all(is.finite(z)))
}
write_r <- function(x, path) terra::writeRaster(x, path, overwrite=TRUE, wopt=WOPT)
plot_r <- function(x, path, title) {
  grDevices::png(path, width=1600, height=1300, res=160)
  on.exit(grDevices::dev.off(), add=TRUE)
  terra::plot(x, main=title)
}
find_tif <- function(folder, stem, required=TRUE) {
  fs <- list.files(folder, pattern="\\.tiff?$", full.names=TRUE, ignore.case=TRUE)
  hits <- fs[tolower(tools::file_path_sans_ext(basename(fs))) == tolower(stem)]
  if (length(hits)>1) stop("Ambiguous raster: ", stem, " in ",folder)
  if (!length(hits)) {
    if (required) stop("Missing raster: ",file.path(folder,paste0(stem,".tif")))
    return(NA_character_)
  }
  hits
}
bio_paths <- function(folder) {
  fs <- list.files(file.path(folder,"bio_clim_var"), "\\.tiff?$", full.names=TRUE,
                   ignore.case=TRUE)
  ids <- vapply(basename(fs), function(s) {
    m <- regexec("bio[_-]?([0-9]{1,2})([^0-9]|$)",tolower(s)); z <- regmatches(tolower(s),m)[[1]]
    if (length(z)) as.integer(z[2]) else NA_integer_
  }, integer(1))
  if (anyNA(ids) || anyDuplicated(ids) || !setequal(ids,1:19))
    stop("bio_clim_var must contain exactly 19 single-band BIO1–BIO19 TIFFs; inspect filenames: ",folder)
  setNames(fs[order(ids)],paste0("bio",1:19))
}
extract_raw <- function(r, xy) {
  # Explicit CRS transformation; ID=FALSE avoids the first-predictor deletion bug.
  v <- terra::vect(xy[,c("Longitude","Latitude"),drop=FALSE],
                    geom=c("Longitude","Latitude"), crs="EPSG:4326")
  v <- terra::project(v,terra::crs(r))
  z <- as.data.frame(terra::extract(r,v,ID=FALSE,raw=TRUE))
  names(z) <- names(r)
  z
}
cells_xy <- function(r, xy) {
  v <- terra::vect(xy[,c("Longitude","Latitude"),drop=FALSE],
                    geom=c("Longitude","Latitude"),crs="EPSG:4326")
  terra::cellFromXY(r,terra::crds(terra::project(v,terra::crs(r))))
}
to_lonlat <- function(r, xy) {
  v <- terra::vect(as.matrix(xy),crs=terra::crs(r))
  z <- as.data.frame(terra::crds(terra::project(v,"EPSG:4326")))
  names(z) <- c("Longitude","Latitude"); z
}
sample_xy <- function(r,n,seed) {
  set.seed(seed)
  z <- terra::spatSample(r,size=n,method="random",na.rm=TRUE,xy=TRUE,
                         as.df=TRUE,values=TRUE)
  if (!nrow(z)) stop("No valid cells to sample")
  to_lonlat(r,z[,c("x","y"),drop=FALSE])
}

# 2. BOUNDARIES, RASTER AUDIT AND ALIGNMENT -------------------------------------
load_africa <- function() {
  if (nzchar(CFG$africa_boundary)) {
    a <- sf::st_read(CFG$africa_boundary,quiet=TRUE)
  } else {
    a <- rnaturalearth::ne_countries(scale=110,continent="Africa",returnclass="sf")
  }
  sf::st_transform(sf::st_make_valid(a),4326)
}
read_one <- function(path,nm) {
  r <- terra::rast(path)
  if (terra::nlyr(r)!=1) stop("Expected one band, found ",terra::nlyr(r),": ",path,
    ". Create a documented temporal summary first; bands are not averaged automatically.")
  if (!nzchar(terra::crs(r))) stop("Missing CRS: ",path)
  # Keep the original integer codes, not R's positional factor codes. LULC is
  # encoded explicitly by encode(), so category metadata is removed on this copy.
  if (is.factor(r)) levels(r) <- NULL
  if (nm %in% names(INPUT_MULTIPLIER)) r <- r*INPUT_MULTIPLIER[[nm]]
  if (nm %in% names(INPUT_OFFSET)) r <- r+INPUT_OFFSET[[nm]]
  names(r) <- nm; r
}
raster_audit <- function(paths) {
  do.call(rbind,lapply(names(paths),function(nm) {
    r <- read_one(paths[[nm]],nm)
    z <- as.numeric(unlist(terra::spatSample(r,10000,method="random",na.rm=TRUE)))
    z <- z[is.finite(z)]
    data.frame(variable=nm,file=paths[[nm]],bands=terra::nlyr(r),
      nrow=terra::nrow(r),ncol=terra::ncol(r),res_x=terra::res(r)[1],res_y=terra::res(r)[2],
      lonlat=terra::is.lonlat(r),crs=terra::crs(r,proj=TRUE),
      sample_min=if(length(z)) min(z) else NA,
      sample_max=if(length(z)) max(z) else NA,
      sample_unique=length(unique(z)),stringsAsFactors=FALSE)
  }))
}
align_one <- function(r,template,type,path) {
  # Average continuous fine-grid inputs; modal categories; nearest binary water.
  # Bilinear is for continuous similar-resolution/upsampling operations.
  same_crs <- terra::same.crs(r,template)
  down <- same_crs && all(terra::res(template)>terra::res(r)*1.5)
  method <- if (type=="categorical") {
    if (down) "mode" else "near"
  } else if (type=="binary") "near" else if (down) "average" else "bilinear"
  if (terra::compareGeom(r,template,stopOnError=FALSE)) {
    r <- terra::mask(r,template)
    return(write_r(r,path))
  }
  # project(template) sets CRS, resolution, origin and extent in one operation.
  rr <- terra::project(r,template,method=method,filename=path,overwrite=TRUE,wopt=WOPT)
  rr
}
build_environment <- function(folder,africa,dir,template=NULL,current=NULL) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  bp <- bio_paths(folder)
  if (is.null(template)) {
    ref <- read_one(bp[[1]],"bio1")
    av <- terra::project(terra::vect(africa),terra::crs(ref))
    template <- terra::mask(terra::crop(ref,av),av)
    template <- terra::ifel(!is.na(template),1,NA)
    template <- write_r(template,file.path(dir,"reference_grid.tif"))
  }
  paths <- bp
  for(i in seq_len(nrow(EXTRA))) {
    f <- find_tif(folder,EXTRA$stem[i],required=is.null(current))
    if (!is.na(f)) paths[EXTRA$name[i]] <- f
  }
  csv(raster_audit(paths),file.path(dir,"input_raster_audit.csv"))
  files <- character(); static <- character()
  for(nm in c(names(bp),EXTRA$name)) {
    typ <- if (nm %in% names(bp)) "continuous" else EXTRA$type[match(nm,EXTRA$name)]
    if (!nm %in% names(paths)) {
      if (!CFG$allow_static_future_nonbioclim) stop("Missing future predictor ",nm,
        ". Supply it or explicitly enable static non-bioclimatic predictors.")
      use <- if(nm=="aspect") c("northness","eastness") else nm
      for(u in use) {
        dest <- file.path(dir,paste0(u,".tif")); write_r(current[[u]],dest); files[u] <- dest
      }
      static <- c(static,nm); next
    }
    r <- read_one(paths[[nm]],nm)
    if (nm=="water") {
      if (CFG$water_nonwater_is_na) r <- terra::ifel(is.na(r),0,r)
      zz <- as.numeric(unlist(terra::spatSample(r,20000,method="random",na.rm=TRUE)))
      if (any(!zz %in% c(0,1))) stop("WaterRasterized must be 0/1; recode documented source classes first.")
      if (!CFG$water_nonwater_is_na && length(unique(zz))==1 && unique(zz)==1)
        stop("WaterRasterized appears to have only water=1 and NA elsewhere. Confirm whether NA means land; then set water_nonwater_is_na=TRUE if appropriate.")
    }
    if (nm=="aspect") {
      angle <- if(CFG$aspect_units=="degrees") r*pi/180 else r
      if (!CFG$aspect_units %in% c("degrees","radians")) stop("Invalid aspect_units")
      for(u in c("northness","eastness")) {
        val <- if(u=="northness") cos(angle) else sin(angle)
        val <- terra::ifel(r==CFG$aspect_flat_code,0,val)
        names(val) <- u; dest <- file.path(dir,paste0(u,".tif"))
        align_one(val,template,"continuous",dest); files[u] <- dest
      }
    } else {
      dest <- file.path(dir,paste0(nm,".tif"))
      align_one(r,template,typ,dest); files[nm] <- dest
    }
  }
  env <- terra::rast(unname(files)); names(env) <- names(files)
  # Flat cells have no preferred exposure; preserve NA where slope is unknown.
  for(u in c("northness","eastness")) {
    z <- terra::ifel(env[["slope"]]==0,0,env[[u]])
    newpath <- file.path(dir,paste0(u,"_flat_corrected.tif")); write_r(z,newpath)
    files[u] <- newpath
  }
  env <- terra::rast(unname(files)); names(env) <- names(files)
  csv(data.frame(variable=names(files),aligned_file=unname(files)),file.path(dir,"aligned_manifest.csv"))
  csv(data.frame(held_constant=static),file.path(dir,"static_future_assumptions.csv"))
  # Uniform coverage across all candidates. This can remove cells; report it.
  valid <- terra::app(env,fun=function(x) {
    if(is.null(dim(x))) return(if(all(is.finite(x))) 1 else NA_real_)
    ifelse(rowSums(!is.finite(x))==0,1,NA_real_)
  },filename=file.path(dir,"complete_predictor_mask.tif"),overwrite=TRUE,wopt=WOPT)
  valid <- terra::mask(valid,template)
  valid <- write_r(valid,file.path(dir,"valid_africa_mask.tif"))
  plot_r(valid,file.path(dir,"valid_predictor_coverage.png"),"Complete candidate-predictor coverage")
  list(env=env,template=template,valid=valid)
}

# 3. OCCURRENCE AUDIT AND SPECIES LABELS ----------------------------------------
clean_occurrences <- function(sp,file,env,africa,dir) {
  d <- read.csv(file,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  needed <- c("Latitude","Longitude","Species")
  if(!all(needed %in% names(d))) stop("Missing columns in ",file)
  # Ignore the empty extra columns in Arabiensis.csv; preserve core source fields.
  d <- d[,needed,drop=FALSE]; d$source_row <- seq_len(nrow(d))
  d$original_label <- as.character(d$Species)
  label <- tolower(trimws(sub("^anopheles[[:space:]]+","",tolower(d$Species))))
  d$exclusion <- ""
  d$exclusion[!label %in% ACCEPTED_LABELS[[sp]]] <- "taxonomic_label_requires_review"
  d$Latitude <- suppressWarnings(as.numeric(d$Latitude))
  d$Longitude <- suppressWarnings(as.numeric(d$Longitude))
  good <- is.finite(d$Latitude)&is.finite(d$Longitude)&abs(d$Latitude)<=90&abs(d$Longitude)<=180
  d$exclusion[!good] <- "invalid_coordinates"
  d$exclusion[good & d$Latitude==0 & d$Longitude==0] <- "zero_zero_coordinate_review"
  idx <- which(d$exclusion=="")
  if(length(idx)) {
    pts <- sf::st_as_sf(d[idx,],coords=c("Longitude","Latitude"),crs=4326)
    inside <- lengths(sf::st_intersects(pts,africa))>0
    d$exclusion[idx[!inside]] <- "outside_Africa_boundary_review"
  }
  idx <- which(d$exclusion=="")
  dup <- duplicated(d[idx,c("Longitude","Latitude")])
  d$exclusion[idx[dup]] <- "exact_coordinate_duplicate"
  idx <- which(d$exclusion=="")
  if(length(idx)) {
    vals <- extract_raw(env,d[idx,]); ok <- finite_rows(vals)
    d$exclusion[idx[!ok]] <- "missing_candidate_predictor"
  }
  idx <- which(d$exclusion=="")
  d$cell <- NA_real_; d$cell[idx] <- cells_xy(env[[1]],d[idx,])
  # Deterministic choice: first source row in each cell; original table is retained.
  dup <- duplicated(d$cell[idx]); d$exclusion[idx[dup]] <- "same_modelling_cell"
  d$Species <- paste("Anopheles",sp)
  csv(d,file.path(dir,"occurrence_audit_all_rows.csv"))
  csv(d[d$exclusion!="",],file.path(dir,"excluded_occurrences_review.csv"))
  z <- d[d$exclusion=="",]; csv(z,file.path(dir,"clean_occurrences.csv"))
  csv(as.data.frame(table(d$original_label)),file.path(dir,"original_label_counts.csv"))
  if(!nrow(z)) stop("No usable occurrences for ",sp)
  z
}

# 4. BACKGROUND AND SPATIAL FOLDS ----------------------------------------------
domain_mask <- function(sp,valid) {
  if(CFG$calibration_mode=="africa") return(valid)
  if(CFG$calibration_mode!="polygons") stop("Invalid calibration_mode")
  path <- CALIBRATION_FILES[[sp]]
  if(!nzchar(path)||!file.exists(path)) stop("Missing calibration polygon for ",sp)
  p <- terra::project(terra::vect(path),terra::crs(valid))
  terra::mask(valid,p)
}
make_background <- function(sp,pres,all_occ,valid,env,seed) {
  dm <- domain_mask(sp,valid)
  if(CFG$background_method=="random") {
    z <- sample_xy(dm,CFG$background_n+nrow(pres)+1000L,seed)
  } else if(CFG$background_method=="target_group") {
    if(!CFG$target_group_effort_confirmed) stop("Target-group background requires confirmed comparable survey effort.")
    z <- do.call(rbind,all_occ[setdiff(names(all_occ),sp)])
    z <- z[,c("Longitude","Latitude"),drop=FALSE]
    z <- z[is.finite(extract_raw(dm,z)[,1]),,drop=FALSE]
  } else stop("Unknown background_method")
  z$cell <- cells_xy(env[[1]],z)
  z <- z[!duplicated(z$cell)&!z$cell %in% pres$cell,,drop=FALSE]
  if(nrow(z)>CFG$background_n) {set.seed(seed);z <- z[sample(seq_len(nrow(z)),CFG$background_n),]}
  if(nrow(z)<CFG$min_test_background*CFG$outer_k) stop("Too few available background cells for ",sp)
  z
}
add_blocks <- function(d) {
  p <- sf::st_as_sf(d,coords=c("Longitude","Latitude"),crs=4326)
  xy <- sf::st_coordinates(sf::st_transform(p,AEA))
  d$block <- paste(floor(xy[,1]/(CFG$block_km*1000)),
                   floor(xy[,2]/(CFG$block_km*1000)),sep="_")
  d
}
spatial_folds <- function(d,k,seed) {
  # A whole block, with BOTH classes, goes to one fold. No nearest-presence mapping.
  ids <- unique(d$block); ix <- match(d$block,ids)
  counts <- cbind(p=tabulate(ix[d$pa==1],length(ids)),b=tabulate(ix[d$pa==0],length(ids)))
  if(sum(counts[,1]>0)<k) stop("Too few occupied blocks for ",k," folds; inspect block size / species coverage")
  set.seed(seed); best <- NULL; best_score <- Inf
  for(i in seq_len(300)) {
    a <- sample(rep(seq_len(k),length.out=length(ids)))
    f <- a[ix]
    ct <- table(factor(f,levels=seq_len(k)),factor(d$pa,levels=0:1))
    if(any(ct[,2]<CFG$min_test_presences)||any(ct[,1]<CFG$min_test_background)) next
    target <- matrix(rep(colSums(ct)/k,each=k),nrow=k)
    score <- sum((ct-target)^2/pmax(target,1)^2)
    if(score<best_score) {best <- f;best_score <- score}
  }
  if(is.null(best)) stop("Could not form balanced spatial folds. Inspect counts; change block size or k explicitly.")
  best
}
save_fold_map <- function(d,path,title) {
  grDevices::png(path,width=1500,height=1300,res=150);on.exit(grDevices::dev.off())
  plot(d$Longitude,d$Latitude,col=d$fold,pch=ifelse(d$pa==1,19,1),cex=0.45,
       xlab="Longitude",ylab="Latitude",main=title)
  legend("bottomleft",legend=paste("Fold",sort(unique(d$fold))),
         col=sort(unique(d$fold)),pch=19,bty="n")
}

# 5. TRAINING-ONLY VARIABLE SELECTION AND ENCODING ------------------------------
select_variables <- function(x) {
  categorical <- intersect("lulc",names(x))
  num <- setdiff(names(x),categorical)
  num <- num[vapply(x[num],function(z) length(unique(z[is.finite(z)]))>1,logical(1))]
  if(length(num)>1) {
    repeat {
      cm <- abs(cor(x[,num,drop=FALSE]));diag(cm) <- 0
      if(max(cm)<=CFG$correlation_limit || length(num)<2) break
      ij <- which(cm==max(cm),arr.ind=TRUE)[1,]
      # Remove the more globally redundant member; deterministic tie break.
      avg <- rowMeans(cm);drop <- if(avg[ij[1]]>=avg[ij[2]]) ij[1] else ij[2]
      num <- num[-drop]
      if(length(num)<2) break
    }
  }
  calc_vif <- function(vars) {
    if(length(vars)<2) return(setNames(rep(1,length(vars)),vars))
    setNames(vapply(seq_along(vars),function(i) {
      yy <- x[[vars[i]]];xx <- as.matrix(x[,vars[-i],drop=FALSE])
      fit <- stats::lm.fit(cbind(1,xx),yy)
      r2 <- 1-sum(fit$residuals^2)/sum((yy-mean(yy))^2)
      if(!is.finite(r2)||r2>=1-1e-10) Inf else 1/(1-r2)
    },numeric(1)),vars)
  }
  repeat {
    vv <- calc_vif(num)
    if(length(num)<2||max(vv)<=CFG$vif_limit) break
    num <- setdiff(num,names(which.max(vv)))
  }
  if(length(categorical) && length(unique(x[[categorical]]))<2) categorical <- character()
  if(!length(c(num,categorical))) stop("No variable predictors remain")
  list(vars=c(num,categorical),vif=calc_vif(num))
}
make_recipe <- function(x,seed) {
  set.seed(seed)
  idx <- if(nrow(x)>CFG$selection_sample_n) sample(seq_len(nrow(x)),CFG$selection_sample_n) else seq_len(nrow(x))
  sel <- select_variables(x[idx,,drop=FALSE])
  lev <- if("lulc" %in% sel$vars) sort(unique(as.numeric(x$lulc))) else numeric()
  list(vars=sel$vars,vif=sel$vif,landcover_levels=lev)
}
encode <- function(x,recipe) {
  z <- x[,setdiff(recipe$vars,"lulc"),drop=FALSE]
  if("lulc" %in% recipe$vars) {
    v <- as.numeric(x$lulc);lev <- recipe$landcover_levels
    # Treatment coding; reference is the first observed training category.
    for(l in lev[-1]) z[[paste0("lc_",safe_id(l))]] <- as.numeric(v==l)
    unknown <- is.na(v)|!v %in% lev
    if(any(unknown)) z[unknown,] <- NA_real_
  }
  as.data.frame(z)
}

# 6. MODEL FITTING AND PREDICTION ----------------------------------------------
grids <- list(
  MAXNET=expand.grid(classes=c("lq","lqh"),regmult=c(1,2,4),stringsAsFactors=FALSE),
  RF=expand.grid(mtry_fraction=c(0.33,0.67),min_node=c(5,15)),
  XGB=expand.grid(depth=c(2,4),eta=c(0.03,0.1),rounds=c(150,350)),
  GAM=expand.grid(k=c(4,6),gamma=c(1,1.4))
)
fit_algorithm <- function(x,y,alg,param,seed) {
  set.seed(seed)
  ratio <- sum(y==0)/sum(y==1)
  if(!is.finite(ratio)||ratio<=0) stop("Both response classes required")
  if(alg=="MAXNET") return(maxnet::maxnet(y,x,
      f=maxnet::maxnet.formula(y,x,classes=as.character(param$classes)),regmult=param$regmult))
  if(alg=="RF") return(ranger::ranger(x=x,y=factor(y,levels=c(0,1)),
      probability=TRUE,num.trees=500,mtry=max(1,round(ncol(x)*param$mtry_fraction)),
      min.node.size=param$min_node,class.weights=c(1,ratio),
      num.threads=CFG$threads,seed=seed))
  if(alg=="XGB") return(xgboost::xgb.train(
      params=list(objective="binary:logistic",eval_metric="auc",max_depth=param$depth,
                  eta=param$eta,subsample=0.8,colsample_bytree=0.8,
                  scale_pos_weight=ratio,nthread=CFG$threads,seed=seed),
      data=xgboost::xgb.DMatrix(as.matrix(x),label=y),nrounds=param$rounds,verbose=0))
  if(alg=="GAM") {
    terms <- vapply(names(x),function(v) {
      nu <- length(unique(x[[v]]))
      if(nu<=5||startsWith(v,"lc_")) v else
        sprintf("s(%s,k=%d,bs='ts')",v,min(param$k,nu-1L))
    },character(1))
    form <- stats::as.formula(paste("y ~",paste(terms,collapse=" + ")),
                              env=asNamespace("mgcv"))
    return(mgcv::gam(form,data=cbind(y=y,x),family=stats::binomial(),
      weights=ifelse(y==1,ratio,1),method="REML",select=TRUE,gamma=param$gamma))
  }
  stop("Unknown algorithm")
}
predict_algorithm <- function(model,x,alg) {
  if(alg=="MAXNET") return(as.numeric(predict(model,x,type="cloglog",clamp=TRUE)))
  if(alg=="RF") return(as.numeric(predict(model,data=x,num.threads=CFG$threads)$predictions[,"1"]))
  if(alg=="XGB") return(as.numeric(predict(model,xgboost::xgb.DMatrix(as.matrix(x)))))
  as.numeric(predict(model,newdata=x,type="response"))
}
fit_bundle <- function(raw,y,alg,param,recipe,seed) {
  x <- encode(raw,recipe);ok <- finite_rows(x)
  if(!all(ok)) stop("Non-finite training predictors after encoding")
  list(model=fit_algorithm(x,y,alg,param,seed),recipe=recipe,algorithm=alg,param=param)
}
predict_bundle <- function(bundle,raw) {
  x <- encode(raw,bundle$recipe);ok <- finite_rows(x);p <- rep(NA_real_,nrow(x))
  if(any(ok)) p[ok] <- predict_algorithm(bundle$model,x[ok,,drop=FALSE],bundle$algorithm)
  p[!is.finite(p)] <- NA_real_;pmin(1,pmax(0,p))
}

# 7. METRICS AND INNER SPATIAL TUNING ------------------------------------------
auc_fixed <- function(y,p) {
  ok <- is.finite(p);if(length(unique(y[ok]))<2) return(NA_real_)
  # Fixed orientation prevents automatic reversal of poor predictions.
  as.numeric(pROC::auc(pROC::roc(y[ok],p[ok],levels=c(0,1),direction="<",quiet=TRUE)))
}
threshold_tss <- function(y,p) {
  ok <- is.finite(p);y <- y[ok];p <- p[ok]
  if(length(unique(y))<2) return(NA_real_)
  ord <- order(p,decreasing=TRUE);y <- y[ord];p <- p[ord]
  ends <- c(which(diff(p)!=0),length(p))
  tp <- cumsum(y)[ends];fp <- ends-tp
  score <- tp/sum(y==1)-fp/sum(y==0)
  # Highest threshold on ties, with prediction rule p >= threshold.
  p[ends[which.max(score)]]
}
metrics <- function(y,p,threshold) {
  coverage <- mean(is.finite(p));ok <- is.finite(p);y <- y[ok];p <- p[ok]
  if(length(unique(y))<2||!is.finite(threshold))
    return(data.frame(AUC=NA,TSS=NA,Sensitivity=NA,Specificity=NA,
      Balanced_accuracy=NA,Boyce=NA,Coverage=coverage,N=length(y),Threshold=threshold))
  b <- p>=threshold
  sens <- mean(b[y==1]);spec <- mean(!b[y==0])
  boyce <- tryCatch(as.numeric(ecospat::ecospat.boyce(
      fit=p[y==0],obs=p[y==1],PEplot=FALSE)$Spearman.cor),error=function(e) NA_real_)
  data.frame(AUC=auc_fixed(y,p),TSS=sens+spec-1,Sensitivity=sens,Specificity=spec,
    Balanced_accuracy=(sens+spec)/2,Boyce=boyce,Coverage=coverage,N=length(y),Threshold=threshold)
}
inner_tune <- function(d,predictors,seed,dir) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  folds <- spatial_folds(d,CFG$inner_k,seed)
  csv(data.frame(row=seq_len(nrow(d)),block=d$block,pa=d$pa,fold=folds),
      file.path(dir,"inner_folds.csv"))
  recipes <- lapply(seq_len(CFG$inner_k),function(f)
    make_recipe(d[folds!=f,predictors,drop=FALSE],seed+f))
  selected <- list();oof <- list();tuning <- list();counts <- 0L
  for(alg in CFG$algorithms) {
    grid <- grids[[alg]];best_score <- -Inf;best_pred <- NULL;best_param <- NULL
    for(j in seq_len(nrow(grid))) {
      pp <- rep(NA_real_,nrow(d));failure <- character()
      for(f in seq_len(CFG$inner_k)) {
        tr <- folds!=f;va <- !tr
        pp[va] <- tryCatch({
          mod <- fit_bundle(d[tr,predictors,drop=FALSE],d$pa[tr],alg,
                            grid[j,,drop=FALSE],recipes[[f]],seed+100*j+f)
          predict_bundle(mod,d[va,predictors,drop=FALSE])
        },error=function(e) {failure <<- c(failure,conditionMessage(e));rep(NA_real_,sum(va))})
      }
      fold_auc <- vapply(seq_len(CFG$inner_k),function(f) auc_fixed(d$pa[folds==f],pp[folds==f]),numeric(1))
      coverage <- mean(is.finite(pp))
      score <- if(all(is.finite(fold_auc))&&coverage>=CFG$minimum_inner_coverage) mean(fold_auc) else -Inf
      counts <- counts+1L
      tuning[[counts]] <- data.frame(Algorithm=alg,Candidate=j,Mean_inner_AUC=score,
        Coverage=coverage,Parameters=paste(names(grid),unlist(grid[j,]),sep="=",collapse=";"),
        Errors=paste(unique(failure),collapse=" | "))
      message("  ",alg," candidate ",j,"/",nrow(grid)," inner AUC=",round(score,3))
      if(score>best_score) {best_score<-score;best_pred<-pp;best_param<-grid[j,,drop=FALSE]}
    }
    if(is.null(best_param)) stop("No valid inner-CV candidate for ",alg,
      ". Inspect category coverage, fold size and errors. No silent algorithm substitution.")
    selected[[alg]] <- best_param;oof[[alg]] <- best_pred
  }
  csv(do.call(rbind,tuning),file.path(dir,"hyperparameter_search.csv"))
  pmat <- do.call(cbind,oof);colnames(pmat)<-names(oof)
  common <- apply(pmat,1,function(z) all(is.finite(z)))
  if(mean(common)<CFG$minimum_inner_coverage) stop("Insufficient common inner prediction coverage")
  summary <- do.call(rbind,lapply(colnames(pmat),function(alg) {
    th <- threshold_tss(d$pa[common],pmat[common,alg])
    cbind(Algorithm=alg,metrics(d$pa[common],pmat[common,alg],th))
  }))
  ww <- if(CFG$ensemble_weight=="TSS") pmax(summary$TSS,0) else if(
    CFG$ensemble_weight=="AUC_above_0.5") pmax(summary$AUC-0.5,0) else stop("Unknown weight rule")
  if(any(!is.finite(ww))||sum(ww)<=0) stop("No positive ensemble skill; inspect models instead of forcing equal weights")
  ww <- setNames(ww/sum(ww),summary$Algorithm);summary$Weight <- ww
  ens <- rep(NA_real_,nrow(d));ens[common] <- as.numeric(pmat[common,,drop=FALSE]%*%ww)
  ens_th <- threshold_tss(d$pa[common],ens[common])
  csv(summary,file.path(dir,"training_only_weights_and_thresholds.csv"))
  csv(cbind(d[,c("Longitude","Latitude","pa","block")],fold=folds,pmat,ensemble=ens),
      file.path(dir,"inner_oof_predictions.csv"))
  saveRDS(selected,file.path(dir,"selected_parameters.rds"))
  list(params=selected,weights=ww,algorithm_threshold=setNames(summary$Threshold,summary$Algorithm),
       threshold=ens_th,summary=summary)
}

# 8. OUTER VALIDATION AND FINAL FIT --------------------------------------------
ensemble_predict <- function(models,weights,raw) {
  pm <- vapply(names(models),function(a) predict_bundle(models[[a]],raw),numeric(nrow(raw)))
  pm <- matrix(pm,nrow=nrow(raw),dimnames=list(NULL,names(models)))
  good <- apply(pm,1,function(z) all(is.finite(z)))
  en <- rep(NA_real_,nrow(raw))
  en[good] <- as.numeric(pm[good,,drop=FALSE]%*%weights[colnames(pm)])
  list(individual=pm,ensemble=en)
}
fit_set <- function(d,predictors,tuned,seed) {
  recipe <- make_recipe(d[,predictors,drop=FALSE],seed)
  models <- setNames(lapply(seq_along(CFG$algorithms),function(i) {
    alg <- CFG$algorithms[i]
    fit_bundle(d[,predictors,drop=FALSE],d$pa,alg,tuned$params[[alg]],recipe,seed+i)
  }),CFG$algorithms)
  list(models=models,recipe=recipe)
}
pack_models <- function(models) {
  lapply(models,function(b) {
    if(b$algorithm=="XGB") {
      b$model <- xgboost::xgb.save.raw(b$model,raw_format="ubj");b$raw_xgb <- TRUE
    }
    b
  })
}
unpack_models <- function(models) {
  lapply(models,function(b) {
    if(isTRUE(b$raw_xgb)) {b$model<-xgboost::xgb.load.raw(b$model);b$raw_xgb<-FALSE}
    b
  })
}
permutation_importance <- function(models,weights,test,predictors,seed) {
  raw <- test[,predictors,drop=FALSE];baseline <- ensemble_predict(models,weights,raw)
  base_auc <- auc_fixed(test$pa,baseline$ensemble)
  vars <- unique(unlist(lapply(models,function(b)b$recipe$vars)))
  set.seed(seed);out <- list();ct <- 0L
  # Hold-out importance: permutation decreases in ensemble AUC; negative values
  # are kept. Group northness/eastness to preserve a meaningful joint aspect test.
  groups <- setNames(lapply(vars,identity),vars)
  asp <- intersect(c("northness","eastness"),vars)
  if(length(asp)) {groups[asp]<-NULL;groups$aspect<-asp}
  for(v in names(groups)) for(j in seq_len(CFG$permutation_repeats)) {
    rr <- raw;ord <- sample(seq_len(nrow(rr)))
    rr[,groups[[v]]] <- rr[ord,groups[[v]],drop=FALSE]
    p <- ensemble_predict(models,weights,rr)$ensemble
    common <- is.finite(p)&is.finite(baseline$ensemble)
    delta <- auc_fixed(test$pa[common],baseline$ensemble[common])-auc_fixed(test$pa[common],p[common])
    ct<-ct+1L;out[[ct]]<-data.frame(Variable=v,Repeat=j,AUC_drop=delta,
      Baseline_AUC=base_auc,N_common=sum(common))
  }
  do.call(rbind,out)
}
fit_species <- function(sp,d,predictors,dir,seed) {
  d <- add_blocks(d)
  d$fold <- spatial_folds(d,CFG$outer_k,seed)
  csv(d[,c("Longitude","Latitude","pa","cell","block","fold")],
      file.path(dir,"outer_fold_assignments.csv"))
  csv(as.data.frame(table(Fold=d$fold,Class=d$pa)),file.path(dir,"outer_fold_counts.csv"))
  save_fold_map(d,file.path(dir,"outer_folds.png"),paste(sp,"spatial validation"))
  results<-list();predictions<-list();importance<-list()
  for(f in seq_len(CFG$outer_k)) {
    fd <- file.path(dir,paste0("outer_fold_",f));dir.create(fd,showWarnings=FALSE)
    checkpoint <- file.path(fd,"complete.rds")
    if(file.exists(checkpoint)) {
      done <- readRDS(checkpoint)
      results[[f]]<-done$metrics;predictions[[f]]<-done$predictions;importance[[f]]<-done$importance
      message("Resuming: outer fold ",f," already completed");next
    }
    message("\n",sp,": OUTER FOLD ",f,"/",CFG$outer_k)
    tr <- d[d$fold!=f,,drop=FALSE];te <- d[d$fold==f,,drop=FALSE]
    if(length(intersect(tr$block,te$block))) stop("Spatial block leakage")
    tuned <- inner_tune(tr,predictors,seed+1000*f,file.path(fd,"inner_tuning"))
    fitted <- fit_set(tr,predictors,tuned,seed+2000*f)
    csv(data.frame(Variable=fitted$recipe$vars),file.path(fd,"selected_predictors.csv"))
    pr <- ensemble_predict(fitted$models,tuned$weights,te[,predictors,drop=FALSE])
    common <- is.finite(pr$ensemble)
    if(mean(common)<CFG$minimum_inner_coverage) stop("Low outer prediction coverage in ",sp," fold ",f,
      "; inspect rare/unseen land-cover classes. No silent exclusion of most test points.")
    met <- do.call(rbind,lapply(CFG$algorithms,function(a)
      cbind(Algorithm=a,metrics(te$pa,pr$individual[,a],tuned$algorithm_threshold[[a]]))))
    met <- rbind(met,cbind(Algorithm="ENSEMBLE",metrics(te$pa,pr$ensemble,tuned$threshold)))
    met$Species<-sp;met$Fold<-f
    pp <- cbind(te[,c("Longitude","Latitude","pa","cell","block","fold")],
                pr$individual,ensemble=pr$ensemble,ensemble_threshold=tuned$threshold)
    pp$ensemble_binary <- as.integer(pp$ensemble>=tuned$threshold)
    imp <- permutation_importance(fitted$models,tuned$weights,te,predictors,seed+3000*f)
    imp$Fold<-f;imp$Species<-sp
    csv(met,file.path(fd,"outer_test_metrics.csv"));csv(pp,file.path(fd,"outer_test_predictions.csv"))
    csv(imp,file.path(fd,"heldout_permutation_importance.csv"))
    saveRDS(list(metrics=met,predictions=pp,importance=imp),checkpoint)
    results[[f]]<-met;predictions[[f]]<-pp;importance[[f]]<-imp
    rm(fitted);gc()
  }
  res <- do.call(rbind,results);csv(res,file.path(dir,"spatial_CV_metrics_all_folds.csv"))
  csv(do.call(rbind,predictions),file.path(dir,"outer_OOF_predictions.csv"))
  imp <- do.call(rbind,importance);csv(imp,file.path(dir,"heldout_importance_all_folds.csv"))
  summary <- do.call(rbind,lapply(unique(res$Algorithm),function(a) {
    z<-res[res$Algorithm==a,];cols<-c("AUC","TSS","Sensitivity","Specificity","Balanced_accuracy","Boyce","Coverage")
    data.frame(Species=sp,Algorithm=a,Metric=cols,
      Mean=vapply(z[cols],mean,numeric(1),na.rm=TRUE),
      SD=vapply(z[cols],stats::sd,numeric(1),na.rm=TRUE),
      Folds_available=vapply(z[cols],function(v)sum(is.finite(v)),integer(1)))
  }))
  csv(summary,file.path(dir,"spatial_CV_summary.csv"))
  # Final development uses all data, with new inner spatial tuning. Outer metrics
  # describe the pre-specified procedure, not a second evaluation of this refit.
  final_dir<-file.path(dir,"final_fit");dir.create(final_dir,showWarnings=FALSE)
  checkpoint<-file.path(final_dir,"model_bundle.rds")
  if(file.exists(checkpoint)) {
    final<-readRDS(checkpoint);final$models<-unpack_models(final$models)
  } else {
    tuned<-inner_tune(d,predictors,seed+50000,file.path(final_dir,"inner_tuning"))
    fitted<-fit_set(d,predictors,tuned,seed+60000)
    rng<-lapply(d[,predictors,drop=FALSE],range,na.rm=TRUE)
    final<-list(models=fitted$models,recipe=fitted$recipe,weights=tuned$weights,
      threshold=tuned$threshold,ranges=rng,landcover_levels=sort(unique(d$lulc)),
      species=sp,predictors=predictors)
    saved<-final;saved$models<-pack_models(final$models);saveRDS(saved,checkpoint)
    csv(data.frame(Variable=final$recipe$vars),file.path(final_dir,"selected_predictors.csv"))
    csv(data.frame(Variable=names(final$recipe$vif),VIF=unname(final$recipe$vif)),
        file.path(final_dir,"final_VIF.csv"))
    csv(data.frame(Algorithm=names(final$weights),Weight=final$weights),file.path(final_dir,"ensemble_weights.csv"))
    csv(data.frame(Species=sp,Threshold=final$threshold,Rule="Inner spatial OOF ensemble MaxSensSpec"),
        file.path(final_dir,"binary_threshold.csv"))
    num<-setdiff(final$recipe$vars,"lulc")
    if(length(num)>1) {
      cm<-cor(d[,num,drop=FALSE]);csv(data.frame(Variable=rownames(cm),cm,check.names=FALSE),
                                   file.path(final_dir,"selected_predictor_correlations.csv"))
      grDevices::png(file.path(final_dir,"correlation_heatmap.png"),1800,1800,res=160)
      stats::heatmap(cm,symm=TRUE,scale="none",margins=c(10,10));grDevices::dev.off()
      grDevices::png(file.path(final_dir,"correlation_dendrogram.png"),1800,1100,res=160)
      plot(stats::hclust(stats::as.dist(1-abs(cm))),main="Selected continuous predictors: 1 - |r|")
      grDevices::dev.off()
    }
  }
  list(final=final,data=d,summary=summary)
}

# 9. RESPONSE CURVES AND RASTER PREDICTION --------------------------------------
response_curves <- function(final,d,dir,seed) {
  set.seed(seed);raw<-d[sample(seq_len(nrow(d)),min(300,nrow(d))),final$predictors,drop=FALSE]
  out<-list();ct<-0L
  # Marginal partial dependence, averaged over observed environments. These are
  # descriptive responses, not causal effects; unrealistic combinations can occur.
  for(v in final$recipe$vars) {
    if(v %in% c("northness","eastness")) next # handled jointly as aspect below
    grid<-if(v=="lulc") final$landcover_levels else unique(as.numeric(
      stats::quantile(d[[v]],probs=seq(0.02,0.98,length.out=40),na.rm=TRUE)))
    for(g in grid) {
      rr<-raw;rr[[v]]<-g;pp<-ensemble_predict(final$models,final$weights,rr)
      ct<-ct+1L;out[[ct]]<-data.frame(Variable=v,Value=g,Ensemble=mean(pp$ensemble,na.rm=TRUE))
    }
  }
  if(any(c("northness","eastness") %in% final$recipe$vars)) {
    for(g in seq(0,360,by=10)) {
      rr<-raw;rr$northness<-cos(g*pi/180);rr$eastness<-sin(g*pi/180)
      pp<-ensemble_predict(final$models,final$weights,rr)
      ct<-ct+1L;out[[ct]]<-data.frame(Variable="aspect_degrees",Value=g,Ensemble=mean(pp$ensemble,na.rm=TRUE))
    }
  }
  z<-do.call(rbind,out);csv(z,file.path(dir,"ensemble_partial_dependence.csv"))
  grDevices::pdf(file.path(dir,"ensemble_partial_dependence.pdf"),width=8,height=6)
  on.exit(grDevices::dev.off(),add=TRUE)
  for(v in unique(z$Variable)) {
    zz<-z[z$Variable==v,];plot(zz$Value,zz$Ensemble,type=if(v=="lulc")"h" else "l",
      xlab=v,ylab="Mean relative suitability",main=paste(final$species,v),ylim=c(0,1))
  }
}
map_predict <- function(final,env,valid,dir) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  suitfile<-file.path(dir,"ensemble_suitability.tif");binfile<-file.path(dir,"binary_suitability.tif")
  done<-file.path(dir,"maps_complete.rds")
  if(file.exists(done)) return(readRDS(done))
  # terra processes file-backed rasters in chunks; no full Africa values() call.
  input<-terra::mask(env,valid)
  multi<-terra::predict(input,final,fun=function(model,data,...) {
    data<-as.data.frame(data);names(data)<-model$predictors
    p<-ensemble_predict(model$models,model$weights,data)
    deviation<-sweep(p$individual,1,p$ensemble,"-")
    sd_alg<-sqrt(as.numeric((deviation^2)%*%model$weights[colnames(p$individual)]))
    used<-model$recipe$vars;continuous<-setdiff(used,"lulc")
    novel<-rep(FALSE,nrow(data))
    for(v in continuous) novel<-novel|data[[v]]<model$ranges[[v]][1]|data[[v]]>model$ranges[[v]][2]
    if("lulc" %in% used) novel<-novel|!data$lulc %in% model$landcover_levels
    # Diagnostic is a univariate range flag, NOT MESS/MOP and not a confidence CI.
    cbind(p$individual,ensemble=p$ensemble,algorithm_SD=sd_alg,novel_environment=as.integer(novel))
  },na.rm=TRUE,cores=1,filename=file.path(dir,"algorithm_and_ensemble_predictions.tif"),
    overwrite=TRUE,wopt=WOPT)
  names(multi)<-c(names(final$models),"ensemble","algorithm_SD","novel_environment")
  write_r(multi[["ensemble"]],suitfile)
  write_r(multi[["ensemble"]]>=final$threshold,binfile)
  write_r(multi[["algorithm_SD"]],file.path(dir,"algorithm_disagreement_SD.tif"))
  write_r(multi[["novel_environment"]],file.path(dir,"environmental_novelty_flag.tif"))
  plot_r(multi[["ensemble"]],file.path(dir,"ensemble_suitability.png"),paste(final$species,"relative suitability"))
  plot_r(multi[["algorithm_SD"]],file.path(dir,"algorithm_disagreement.png"),paste(final$species,"algorithm disagreement"))
  b<-terra::rast(binfile);a<-terra::cellSize(b,unit="km")
  csv(data.frame(Species=final$species,Threshold=final$threshold,
      Suitable_km2=terra::global(b*a,"sum",na.rm=TRUE)[1,1],
      Modelled_km2=terra::global(terra::ifel(!is.na(b),a,NA),"sum",na.rm=TRUE)[1,1]),
      file.path(dir,"suitable_area.csv"))
  paths<-list(suitability=suitfile,binary=binfile,dir=dir)
  saveRDS(paths,done);paths
}

# 10. FUTURE PROJECTIONS, GCM AGREEMENT AND RANGE CHANGE -------------------------
shift_map <- function(current_binary,future_binary,path) {
  # 0=stable unsuitable; 1=loss; 2=gain; 3=stable suitable.
  # Encoding b0 + 2*b1 distinguishes both kinds of stability.
  write_r(current_binary+2*future_binary,path)
}
future_species <- function(final,current,scenarios,dir) {
  if(!length(scenarios)) return(invisible(NULL))
  all_maps<-list()
  for(id in names(scenarios)) {
    s<-scenarios[[id]];dd<-file.path(dir,"future",id)
    all_maps[[id]]<-map_predict(final,s$env,s$valid,dd)
    b0<-terra::rast(current$binary);b1<-terra::rast(all_maps[[id]]$binary)
    shift_map(b0,b1,file.path(dd,"range_shift_four_classes.tif"))
  }
  groups<-split(seq_len(nrow(FUTURE)),paste(FUTURE$ssp,FUTURE$period,sep="_"))
  ids<-vapply(seq_len(nrow(FUTURE)),function(i)safe_id(paste(FUTURE$gcm[i],FUTURE$ssp[i],FUTURE$period[i],sep="_")),character(1))
  for(g in names(groups)) {
    ii<-groups[[g]];dd<-file.path(dir,"future_consensus",safe_id(g))
    dir.create(dd,recursive=TRUE,showWarnings=FALSE)
    ss<-terra::rast(vapply(all_maps[ids[ii]],`[[`,character(1),"suitability"))
    bb<-terra::rast(vapply(all_maps[ids[ii]],`[[`,character(1),"binary"))
    # No na.rm=TRUE: consensus requires all specified GCMs at the pixel.
    mn<-terra::app(ss,mean,na.rm=FALSE);frac<-terra::app(bb,mean,na.rm=FALSE)
    spread<-if(length(ii)>1) terra::app(ss,sd,na.rm=FALSE) else mn*NA_real_
    write_r(mn,file.path(dd,"GCM_mean_suitability.tif"))
    write_r(spread,file.path(dd,"GCM_suitability_SD.tif"))
    write_r(frac,file.path(dd,"fraction_GCMs_predicting_suitable.tif"))
    # Agreement regardless of class is separately labelled, never 'robust suitability'.
    write_r(terra::ifel(frac>=0.5,frac,1-frac),file.path(dd,"majority_class_agreement.tif"))
    binary<-write_r(mn>=final$threshold,file.path(dd,"binary_suitability.tif"))
    shifts<-shift_map(terra::rast(current$binary),binary,file.path(dd,"range_shift_four_classes.tif"))
    area<-terra::cellSize(shifts,unit="km")
    stats<-do.call(rbind,lapply(0:3,function(cl) data.frame(Code=cl,
      Class=c("Stable unsuitable","Loss","Gain","Stable suitable")[cl+1],
      Area_km2=terra::global(terra::ifel(shifts==cl,area,0),"sum",na.rm=TRUE)[1,1])))
    csv(stats,file.path(dd,"range_shift_areas.csv"))
    plot_r(shifts,file.path(dd,"range_shift.png"),paste(final$species,g,"0 absent; 1 loss; 2 gain; 3 retained"))
    csv(FUTURE[ii,,drop=FALSE],file.path(dd,"GCM_members.csv"))
  }
}

# 11. GEOGRAPHIC OVERLAP AND POTENTIAL SPECIES RICHNESS --------------------------
geographic_overlap <- function(binary_files,dir) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  if(length(binary_files)!=nrow(SPECIES)) stop("All nine species are needed for the planned richness/overlap outputs")
  r<-terra::rast(unname(binary_files));names(r)<-names(binary_files)
  richness<-terra::app(r,sum,na.rm=FALSE) # any unmodelled species -> NA, not zero
  write_r(richness,file.path(dir,"potential_species_richness_0_to_9.tif"))
  plot_r(richness,file.path(dir,"potential_species_richness.png"),"Potential richness of the nine modelled Anopheles species")
  common_n<-terra::app(r,function(x) {
    if(is.null(dim(x))) return(sum(is.finite(x)))
    rowSums(is.finite(x))
  })
  write_r(common_n,file.path(dir,"number_of_species_with_predictions.tif"))
  pairs<-combn(seq_len(terra::nlyr(r)),2);out<-vector("list",ncol(pairs))
  area<-terra::cellSize(r[[1]],unit="km")
  for(j in seq_len(ncol(pairs))) {
    i<-pairs[1,j];k<-pairs[2,j];a<-r[[i]];b<-r[[k]]
    valid<-!is.na(a)&!is.na(b)
    aa<-terra::global(terra::ifel(valid & a==1,area,NA),"sum",na.rm=TRUE)[1,1]
    ab<-terra::global(terra::ifel(valid & b==1,area,NA),"sum",na.rm=TRUE)[1,1]
    shared<-terra::global(terra::ifel(valid & a==1 & b==1,area,NA),"sum",na.rm=TRUE)[1,1]
    aa<-if(is.finite(aa))aa else 0;ab<-if(is.finite(ab))ab else 0
    shared<-if(is.finite(shared))shared else 0
    union<-aa+ab-shared
    domain_area<-terra::global(terra::ifel(valid,area,NA),"sum",na.rm=TRUE)[1,1]
    out[[j]]<-data.frame(Species_1=names(r)[i],Species_2=names(r)[k],
      Shared_km2=shared,Union_km2=union,Common_domain_km2=domain_area,
      Jaccard=if(union>0)shared/union else NA_real_,
      Sorensen=if(aa+ab>0)2*shared/(aa+ab) else NA_real_)
    # 0=neither, 1=species1 only, 2=species2 only, 3=both; NA preserved.
    write_r(a+2*b,file.path(dir,paste0(names(r)[i],"__",names(r)[k],"_overlap.tif")))
  }
  z<-do.call(rbind,out);csv(z,file.path(dir,"pairwise_geographic_overlap_36_pairs.csv"))
  for(metric in c("Jaccard","Sorensen")) {
    m<-matrix(NA_real_,nrow(SPECIES),nrow(SPECIES),dimnames=list(names(r),names(r)));diag(m)<-1
    for(j in seq_len(nrow(z))) {m[z$Species_1[j],z$Species_2[j]]<-z[[metric]][j];m[z$Species_2[j],z$Species_1[j]]<-z[[metric]][j]}
    csv(data.frame(Species=rownames(m),m,check.names=FALSE),file.path(dir,paste0(metric,"_matrix.csv")))
    grDevices::png(file.path(dir,paste0(metric,"_heatmap.png")),1800,1600,res=160)
    stats::heatmap(m,Rowv=NA,Colv=NA,scale="none",margins=c(12,12),main=paste(metric,"area-weighted geographic overlap"))
    grDevices::dev.off()
  }
}

# 12. COMMON PCA-ENV NICHE OVERLAP ---------------------------------------------
niche_overlap <- function(env,valid,occ,dir,seed) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  xy<-sample_xy(valid,CFG$background_pca_n,seed)
  bx<-extract_raw(env,xy)
  # Do not put arbitrary categorical codes in PCA. Binary water is excluded from
  # this continuous-environment comparison and remains available to the SDMs.
  bx<-bx[,setdiff(names(bx),c("lulc","water")),drop=FALSE]
  sel<-select_variables(bx)$vars
  pc<-stats::prcomp(bx[,sel,drop=FALSE],center=TRUE,scale.=TRUE)
  if(ncol(pc$x)<2) stop("At least two nonredundant PCA dimensions needed")
  scores<-as.data.frame(pc$x[,1:2,drop=FALSE])
  pres_scores<-lapply(occ,function(d) {
    x<-extract_raw(env,d)
    as.data.frame(predict(pc,x[,sel,drop=FALSE])[,1:2,drop=FALSE])
  })
  # A shared grid extent includes environmental sample AND all occurrence scores.
  glob<-rbind(scores,do.call(rbind,pres_scores))
  grids_env<-list()
  for(sp in names(occ)) {
    bgsp<-scores
    if(CFG$calibration_mode=="polygons") {
      sx<-sample_xy(domain_mask(sp,valid),CFG$background_pca_n,seed+match(sp,names(occ)))
      ex<-extract_raw(env,sx)
      bgsp<-as.data.frame(predict(pc,ex[,sel,drop=FALSE])[,1:2,drop=FALSE])
    }
    grids_env[[sp]]<-ecospat::ecospat.grid.clim.dyn(glob=glob,glob1=bgsp,
      sp=pres_scores[[sp]],R=CFG$niche_grid_resolution,th.sp=0,th.env=0)
    csv(pres_scores[[sp]],file.path(dir,paste0(sp,"_occurrence_PC_scores.csv")))
  }
  pairs<-combn(names(occ),2);out<-vector("list",ncol(pairs))
  for(j in seq_len(ncol(pairs))) {
    ov<-ecospat::ecospat.niche.overlap(grids_env[[pairs[1,j]]],grids_env[[pairs[2,j]]],cor=TRUE)
    out[[j]]<-data.frame(Species_1=pairs[1,j],Species_2=pairs[2,j],Schoener_D=ov$D,Hellinger_I=ov$I)
  }
  csv(do.call(rbind,out),file.path(dir,"PCA_env_overlap_36_pairs.csv"))
  csv(data.frame(Variable=rownames(pc$rotation),pc$rotation,check.names=FALSE),file.path(dir,"PCA_loadings.csv"))
  explained<-pc$sdev^2/sum(pc$sdev^2)
  csv(data.frame(Axis=seq_along(explained),Explained_variance=explained,Cumulative=cumsum(explained)),
      file.path(dir,"PCA_explained_variance.csv"))
  saveRDS(list(pca=pc,variables=sel,density_grids=grids_env),file.path(dir,"PCA_env_objects.rds"))
  grDevices::pdf(file.path(dir,"shared_PCA_environment.pdf"),9,7)
  plot(scores,pch=16,col=grDevices::adjustcolor("grey50",alpha.f=0.12),
       xlab=paste0("PC1 (",round(explained[1]*100,1),"%)"),
       ylab=paste0("PC2 (",round(explained[2]*100,1),"%)"),main="Shared continuous environmental space")
  colors<-grDevices::hcl.colors(length(occ),"Dark 3")
  for(i in seq_along(occ)) points(pres_scores[[i]],pch=16,cex=0.35,col=grDevices::adjustcolor(colors[i],0.5))
  legend("topright",legend=names(occ),col=colors,pch=16,cex=0.7,bty="n")
  grDevices::dev.off()
  writeLines(c("Descriptive realized environmental overlap using two common PCA axes.",
    "D and I are corrected for environmental availability within the specified calibration domains.",
    "This correction does not remove all occurrence sampling bias.",
    "No equivalency/similarity randomization tests or future occupied-niche claims are made.",
    paste("Variance represented by PC1+PC2:",round(sum(explained[1:2]),4))),
    file.path(dir,"interpretation.txt"))
}

# 13. RUN FINGERPRINT AND MAIN WORKFLOW ----------------------------------------
check_run <- function() {
  folders<-unique(c(ENV_DIR,OCC_DIR,FUTURE$folder))
  paths<-unique(unlist(lapply(folders,function(p)list.files(p,recursive=TRUE,full.names=TRUE,
                                        pattern="\\.(tif|tiff|csv)$",ignore.case=TRUE))))
  extra_paths<-c(CFG$africa_boundary,unname(CALIBRATION_FILES))
  # Include all shapefile companion files, if applicable.
  for(p in extra_paths[nzchar(extra_paths)]) {
    if(grepl("\\.shp$",p,ignore.case=TRUE)) {
      stem<-tools::file_path_sans_ext(basename(p))
      fs<-list.files(dirname(p),full.names=TRUE)
      paths<-c(paths,fs[tools::file_path_sans_ext(basename(fs))==stem])
    } else paths<-c(paths,p)
  }
  paths<-sort(unique(paths));info<-file.info(paths)
  manifest<-data.frame(path=paths,size=info$size,mtime=as.character(info$mtime),stringsAsFactors=FALSE)
  # Fingerprint function bodies too, so editing implementation invalidates caches.
  funcs<-ls(envir=.GlobalEnv)
  funcs<-funcs[vapply(funcs,function(n)is.function(get(n,envir=.GlobalEnv)),logical(1))]
  code<-lapply(funcs,function(n)paste(deparse(body(get(n,envir=.GlobalEnv))),collapse="\n"))
  versions<-setNames(vapply(required,function(p)as.character(utils::packageVersion(p)),character(1)),required)
  fingerprint<-digest::digest(list(CFG,SPECIES,ACCEPTED_LABELS,EXTRA,INPUT_MULTIPLIER,
    INPUT_OFFSET,FUTURE,CALIBRATION_FILES,manifest,code,versions),algo="sha256")
  path<-file.path(OUT_DIR,"run_fingerprint.rds")
  if(file.exists(path)&&!identical(readRDS(path),fingerprint))
    stop("Inputs/settings/code/package versions changed. Choose a NEW OUT_DIR to prevent stale checkpoints.")
  saveRDS(fingerprint,path);csv(manifest,file.path(OUT_DIR,"input_file_manifest.csv"))
  saveRDS(list(CFG=CFG,SPECIES=SPECIES,EXTRA=EXTRA,FUTURE=FUTURE,
    ACCEPTED_LABELS=ACCEPTED_LABELS,CALIBRATION_FILES=CALIBRATION_FILES),file.path(OUT_DIR,"configuration.rds"))
  writeLines(capture.output(sessionInfo()),file.path(OUT_DIR,"sessionInfo.txt"))
  fingerprint
}
main <- function() {
  if(!dir.exists(ENV_DIR)||!dir.exists(OCC_DIR)) stop("Edit ROOT / ENV_DIR / OCC_DIR first")
  if(!all(file.exists(file.path(OCC_DIR,SPECIES$file)))) stop("Some of the nine occurrence CSVs are missing")
  if(nrow(FUTURE)&&anyDuplicated(FUTURE[,c("gcm","ssp","period")])) stop("Duplicate future GCM/SSP/period rows")
  if(CFG$background_method=="target_group"&&!CFG$target_group_effort_confirmed)
    stop("Confirm shared sampling effort before using target-group background")
  check_run()
  africa<-load_africa()
  ed<-file.path(OUT_DIR,"environment_current");manifest<-file.path(ed,"aligned_manifest.csv")
  ep<-file.path(ed,"environment_complete.rds")
  if(file.exists(ep)) {
    paths<-read.csv(manifest,stringsAsFactors=FALSE)
    current<-list(env=terra::rast(paths$aligned_file),
      template=terra::rast(file.path(ed,"reference_grid.tif")),
      valid=terra::rast(file.path(ed,"valid_africa_mask.tif")))
    names(current$env)<-paths$variable
  } else {
    current<-build_environment(ENV_DIR,africa,ed)
    saveRDS(TRUE,ep)
  }
  # Re-run this inexpensive audit on every invocation; it never changes originals.
  occ<-list();summary<-list()
  for(i in seq_len(nrow(SPECIES))) {
    sp<-SPECIES$species[i];dd<-file.path(OUT_DIR,sp)
    dir.create(dd,recursive=TRUE,showWarnings=FALSE)
    occ[[sp]]<-clean_occurrences(sp,file.path(OCC_DIR,SPECIES$file[i]),current$env,africa,dd)
    summary[[i]]<-data.frame(Species=sp,Raw_rows=nrow(read.csv(file.path(OCC_DIR,SPECIES$file[i]))),
                            Retained_unique_cells=nrow(occ[[sp]]))
  }
  csv(do.call(rbind,summary),file.path(OUT_DIR,"occurrence_summary.csv"))
  print(do.call(rbind,summary))
  if(RUN_MODE=="prepare") {
    message("\nPREPARATION COMPLETE. Inspect input_raster_audit.csv, coverage PNG and excluded_occurrences_review.csv.",
      "\nTo fit models locally, set RUN_MODE='all' and source this script again.")
    return(invisible(NULL))
  }
  scenarios<-list()
  if(nrow(FUTURE)) for(i in seq_len(nrow(FUTURE))) {
    id<-safe_id(paste(FUTURE$gcm[i],FUTURE$ssp[i],FUTURE$period[i],sep="_"))
    if(id %in% names(scenarios)) stop("Scenario identifier collision")
    fd<-file.path(OUT_DIR,"environment_future",id)
    fm<-file.path(fd,"aligned_manifest.csv");flag<-file.path(fd,"environment_complete.rds")
    if(file.exists(flag)) {
      files<-read.csv(fm,stringsAsFactors=FALSE);rr<-terra::rast(files$aligned_file);names(rr)<-files$variable
      scenarios[[id]]<-list(env=rr,valid=terra::rast(file.path(fd,"valid_africa_mask.tif")))
    } else {
      scenarios[[id]]<-build_environment(FUTURE$folder[i],africa,fd,
        template=current$template,current=current$env)
      saveRDS(TRUE,flag)
    }
  }
  maps<-list();failures<-list();evaluations<-list()
  for(i in seq_len(nrow(SPECIES))) {
    sp<-SPECIES$species[i];dd<-file.path(OUT_DIR,sp);seed<-CFG$seed+100000*i
    message("\n================ ",sp," ================")
    tryCatch({
      pres<-occ[[sp]];dm<-domain_mask(sp,current$valid)
      in_domain<-is.finite(extract_raw(dm,pres)[,1])
      csv(pres[!in_domain,,drop=FALSE],file.path(dd,"presences_outside_calibration_area.csv"))
      pres<-pres[in_domain,,drop=FALSE]
      if(nrow(pres)<CFG$min_species_presences) stop("Too few species-level unique cells after cleaning: ",nrow(pres))
      bg<-make_background(sp,pres,occ,current$valid,current$env,seed)
      csv(bg,file.path(dd,"background_points.csv"))
      xy<-rbind(pres[,c("Longitude","Latitude","cell")],bg[,c("Longitude","Latitude","cell")])
      x<-extract_raw(current$env,xy)
      if(!all(finite_rows(x))) stop("Unexpected missing predictors after background/occurrence filtering")
      d<-cbind(xy,pa=c(rep(1L,nrow(pres)),rep(0L,nrow(bg))),x)
      fit<-fit_species(sp,d,names(current$env),dd,seed)
      evaluations[[sp]]<-fit$summary
      if(CFG$make_response_curves) response_curves(fit$final,fit$data,dd,seed)
      maps[[sp]]<-map_predict(fit$final,current$env,current$valid,file.path(dd,"current_maps"))
      future_species(fit$final,maps[[sp]],scenarios,dd)
      saveRDS(list(species=sp,completed=Sys.time()),file.path(dd,"species_complete.rds"))
      if(file.exists(file.path(dd,"FAILED.txt"))) unlink(file.path(dd,"FAILED.txt"))
      rm(fit);gc()
    },error=function(e) {
      failures[[sp]]<<-data.frame(Species=sp,Error=conditionMessage(e))
      writeLines(conditionMessage(e),file.path(dd,"FAILED.txt"))
      message("FAILED ",sp,": ",conditionMessage(e),"\nContinuing with the other species.")
    })
  }
  if(length(evaluations)) csv(do.call(rbind,evaluations),file.path(OUT_DIR,"all_species_CV_summary.csv"))
  if(length(failures)) {
    csv(do.call(rbind,failures),file.path(OUT_DIR,"failed_species.csv"))
    stop("Some species/stages failed. Inspect failed_species.csv. No incomplete nine-species overlap is reported.")
  }
  if(file.exists(file.path(OUT_DIR,"failed_species.csv"))) unlink(file.path(OUT_DIR,"failed_species.csv"))
  geographic_overlap(setNames(vapply(maps,`[[`,character(1),"binary"),names(maps)),
                      file.path(OUT_DIR,"overlap_current"))
  if(length(scenarios)) for(id in names(scenarios)) {
    bf<-setNames(file.path(OUT_DIR,SPECIES$species,"future",id,"binary_suitability.tif"),SPECIES$species)
    geographic_overlap(bf,file.path(OUT_DIR,"overlap_future",id))
  }
  if(CFG$make_niche_overlap) niche_overlap(current$env,current$valid,occ,
                                          file.path(OUT_DIR,"niche_overlap"),CFG$seed)
  message("\nCOMPLETE. Results: ",normalizePath(OUT_DIR,winslash="/"))
}

# Run ONLY when you source this file on your computer.
main()
