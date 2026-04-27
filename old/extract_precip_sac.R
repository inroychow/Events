setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
# Load packages
pacman::p_load(sf, terra, exactextractr, data.table, dplyr, lubridate, stringr, prism, fixest, future, future.apply, progressr, tidyverse, fst, patchwork, scales)

#Load in base data
htdt = fread("data/htdt_allconus.xls")

nfip = fread("L:\\Wetland Flood Mitigation\\NFIP\\Claims\\FimaNfipClaimsV2.csv")
setDT(nfip)  # raw NFIP table

nfip_panel <- {
  nfip[, `:=`(
    date     = as.Date(dateOfLoss),
    tract_id = as.character(censusTract),
    loss_amt = rowSums(.SD, na.rm = TRUE)
  ), .SDcols = c("netBuildingPaymentAmount",
                 "netContentsPaymentAmount")][
                   , .(loss_amt = sum(loss_amt), n_claims = .N), by = .(tract_id, date)
                 ][
                   CJ(tract_id = unique(tract_id),
                      date     = seq(min(date), max(date), by = "day")),
                   on = .(tract_id, date)
                 ][
                   is.na(loss_amt), `:=`(loss_amt = 0, n_claims = 0)
                 ][
                   , any_claim := as.integer(n_claims > 0)
                 ][]
}



#### EXTRACT PRISM #####

# Choose a local folder 
prism_set_dl_dir("prism_archive")  

# Date range to match nfip panel
min_d <- min(nfip_panel$date)
max_d <- max(nfip_panel$date)

# PRISM daily starts 1981-01-01
min_d <- max(min_d, as.Date("1981-01-01"))

# Download daily precipitation #4km resolution
get_prism_dailys(
  type     = "ppt",
  minDate  = min_d,
  maxDate  = max_d,
  keepZip  = FALSE,
  resolution = "4km"
)

# ----------------------------------------------------------------------------------------
# =============================================================================
# PRISM Spatial Subsetting for Sacramento HUC4 (1802)
# Run this BEFORE the main extraction loop in sacramento_flood_pipeline.R
# It crops each daily PRISM .bil to the Sacramento bounding box and saves
# lightweight regional rasters — cuts extraction time
# =============================================================================

# =============================================================================
# 1. SUBSET WATERSHEDS TO SACRAMENTO HUC4
# =============================================================================

watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

# Base R indexing avoids dplyr choking on list 
huc12_sac <- watersheds[str_starts(watersheds$huc12, "1802"), ]

geom <- st_geometry(huc12_sac)

huc12_sac <- st_drop_geometry(huc12_sac) %>%
  mutate(across(where(is.list),
                ~ vapply(., function(x) paste(unlist(x), collapse = "|"), character(1))))

st_geometry(huc12_sac) <- geom
huc12_sac <- st_as_sf(huc12_sac)

# Drop Z/M dimensions from compound CRS (NAD83 + NAVD88 height -> plain NAD83)
huc12_sac <- huc12_sac %>%
  st_zm(drop = TRUE) %>%
  st_transform(4269)   # NAD83 geographic -- matches PRISM native CRS

cat(sprintf("Sacramento HUC12s: %d polygons\n", nrow(huc12_sac)))

# =============================================================================
# 2. BUILD SACRAMENTO BOUNDING BOX FOR CROPPING
# =============================================================================

PRISM_DIR     <- "prism_archive"
SAC_CACHE_DIR <- "prism_sac"

dir.create(SAC_CACHE_DIR, showWarnings = FALSE)

sac_bbox <- huc12_sac %>%
  st_bbox() %>%
  { ext(.$xmin - 0.1, .$xmax + 0.1, .$ymin - 0.1, .$ymax + 0.1) }

cat(sprintf("Sacramento bounding box (buffered):\n"))
cat(sprintf("  xmin=%.3f  xmax=%.3f  ymin=%.3f  ymax=%.3f\n",
            sac_bbox$xmin, sac_bbox$xmax, sac_bbox$ymin, sac_bbox$ymax))

# =============================================================================
# 3. INDEX EXISTING PRISM ARCHIVE
# =============================================================================

prism_set_dl_dir(PRISM_DIR)

# Find all .bil files  — one per dated subfolder
prism_files <- list.files(
  PRISM_DIR,
  pattern    = "\\.bil$",
  recursive  = TRUE,
  full.names = TRUE
)

# Date is in the folder/filename: prism_ppt_us_25m_YYYYMMDD
prism_dates <- str_extract(basename(prism_files), "\\d{8}") %>%
  as.Date(format = "%Y%m%d")

prism_index <- data.table(file = prism_files, date = prism_dates) %>%
  filter(!is.na(date)) %>%
  arrange(date)

cat(sprintf("Found %d PRISM files: %s to %s\n",
            nrow(prism_index),
            min(prism_index$date),
            max(prism_index$date)))
dir.create(SAC_CACHE_DIR, showWarnings = FALSE)

# =============================================================================
# CROP AND SAVE REGIONAL TILES
# =============================================================================
# Output: prism_sac/PRISM_ppt_4kmD2_YYYYMMDD_SAC1802.tif
# GeoTIFF with DEFLATE compression -- ~5-15KB each vs ~500KB+ for national .bil
# Already-cropped files are skipped 

plan(multisession, workers = 15)

handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total files | ETA: :eta",
  width    = 80,
  complete = "="
))

crop_one <- function(i, to_process, sac_bbox, SAC_CACHE_DIR, p) {
  on.exit(p())  # tick progress bar when function exits, even on error
  
  row  <- to_process[i]
  r_in <- tryCatch(
    terra::rast(row$file),
    error = function(e) {
      message(sprintf("\n  Skipping %s -- could not read: %s", row$file, e$message))
      return(NULL)
    }
  )
  if (is.null(r_in)) return(invisible(NULL))
  
  r_crop   <- terra::crop(r_in, sac_bbox)
  date_str <- format(row$date, "%Y%m%d")
  out_file <- file.path(SAC_CACHE_DIR,
                        sprintf("prism_ppt_us_25m_%s_SAC1802.tif", date_str))
  
  terra::writeRaster(r_crop, out_file, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  
  invisible(out_file)
}

# Run with progress reporting
with_progress({
  p <- progressor(steps = nrow(to_process))
  
  future_lapply(
    seq_len(nrow(to_process)),
    crop_one,
    to_process    = to_process,
    sac_bbox      = sac_bbox,
    SAC_CACHE_DIR = SAC_CACHE_DIR,
    p             = p,
    future.seed   = TRUE   
  )
})

plan(sequential)

cat(sprintf("\nDone. Cropped tiles saved to %s/\n", SAC_CACHE_DIR))


##############

# =============================================================================
# EXTRACT HUC12-LEVEL DAILY PRECIP FROM CROPPED PRISM TILES
# =============================================================================
tif_files <- list.files("prism_sac", pattern = "\\.tif$", full.names = TRUE)
dates     <- as.Date(str_extract(basename(tif_files), "\\d{8}"), "%Y%m%d")

tif_dt <- data.table(file = tif_files, date = dates)
tif_dt[, year := year(date)]
years <- sort(unique(tif_dt$year))

PRECIP_HUC12_DIR <- "data/precip_huc12_tiles"
dir.create(PRECIP_HUC12_DIR, showWarnings = FALSE, recursive = TRUE)
tif_dt <- data.table(file = tif_files, date = dates)
tif_dt[, yr_int := data.table::year(date)]
years <- sort(unique(tif_dt$yr_int))
already_done <- list.files(PRECIP_HUC12_DIR, pattern = "\\.fst$")

process_year <- function(yr, tif_dt, huc12_sac, out_dir, already_done) {
  yr <- as.integer(yr)
  yr_rows <- tif_dt[tif_dt$yr_int == yr, ]
  
  huc12_todo_idx <- which(
    !sprintf("precip_%s_%d.fst", huc12_sac$huc12, yr) %in% already_done
  )
  
  if (length(huc12_todo_idx) == 0) return(invisible(yr))
  
  r_yr <- terra::rast(yr_rows$file)
  
  # ----  reproject huc12 to match raster CRS ----
  huc12_sac <- sf::st_transform(huc12_sac, crs = terra::crs(r_yr))
  # ---------------------------------------------------
  
  wide <- exactextractr::exact_extract(
    r_yr,
    huc12_sac[huc12_todo_idx, ],
    fun      = "mean",
    progress = FALSE
  )
  
  for (j in seq_along(huc12_todo_idx)) {
    i        <- huc12_todo_idx[j]
    huc_id   <- huc12_sac$huc12[i]
    out_file <- file.path(out_dir, sprintf("precip_%s_%d.fst", huc_id, yr))
    result   <- data.table::data.table(
      huc12     = huc_id,
      date      = yr_rows$date,
      precip_mm = as.numeric(wide[j, ])
    )
    fst::write_fst(result, out_file, compress = 75)
  }
  
  invisible(yr)
}

plan(multisession, workers = 15)
handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total years | ETA: :eta",
  width    = 80,
  complete = "="
))
with_progress({
  p <- progressor(steps = length(years))
  future_lapply(years, function(yr) {
    on.exit(p())
    process_year(yr, tif_dt, huc12_sac, PRECIP_HUC12_DIR, already_done)
  }, future.seed = TRUE)
})
plan(sequential)
message("Done. Files written to ", PRECIP_HUC12_DIR)

####

files <- list.files(PRECIP_HUC12_DIR, pattern = "\\.fst$", full.names = TRUE)
precip_all <- rbindlist(lapply(files, function(f) {
  tryCatch(
    read_fst(f),
    error = function(e) {
      warning("Skipping: ", f, " — ", conditionMessage(e))
      NULL
    }
  )
}), fill = TRUE)

head(precip_all)
# saveRDS(precip_all, "data/precip_huc12_tiles/precip_huc12s_all.rds")
