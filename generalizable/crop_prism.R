# =============================================================================
# PRISM SPATIAL SUBSETTING + EXTRACTION FOR MULTIPLE HUC4s
# HUC4s: 1806, 1801, 1805, 1804, 1803
# Assumes PRISM daily ppt already downloaded to "prism_archive"
# =============================================================================

library(sf); library(terra); library(exactextractr)
library(data.table); library(dplyr); library(lubridate)
library(stringr); library(prism); library(future); library(future.apply)
library(progressr); library(fst)

# Get which huc4s to extract from 

htdt <- fread("data/htdt_allconus.xls")
htdt_long <- htdt[, .(
  tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ", "))
), by = huc12]

htdt_long[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
htdt_long[, tract := as.character(trimws(tract))]
setnames(htdt_long, "huc12", "upstream_huc12")
# add huc4 from upstream huc12
htdt_long[, huc4 := str_sub(upstream_huc12, 1, 4)]
# keep only rows where tract belongs to that same huc4
htdt_long <- htdt_long %>%
  inner_join(tract_huc4, by = c("tract", "huc4")) %>%
  distinct(upstream_huc12, tract, huc4)


# TARGET HUC4s ----------------------------------
target_huc4s <- sort(unique(htdt_long$huc4))
# -----------------------------------------------

# =============================================================================
# 1. SUBSET WATERSHEDS TO TARGET HUC4s
# =============================================================================
watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

# Filter to any HUC12 whose prefix matches one of the target HUC4s
huc12_target <- watersheds[str_sub(watersheds$huc12, 1, 4) %in% target_huc4s, ]
huc12_target <- huc12_target %>%
  dplyr::select(huc12, geometry) #make smaller for parallelizing
# Flatten any list columns (same fix as original)
geom <- st_geometry(huc12_target)
huc12_target <- st_drop_geometry(huc12_target) %>%
  mutate(across(where(is.list),
                ~ vapply(., function(x) paste(unlist(x), collapse = "|"), character(1))))
st_geometry(huc12_target) <- geom
huc12_target <- st_as_sf(huc12_target)

# Drop Z/M, reproject to NAD83 to match PRISM native CRS
huc12_target <- huc12_target %>%
  st_zm(drop = TRUE) %>%
  st_transform(4269)

cat(sprintf("Target HUC12s: %d polygons across HUC4s: %s\n",
            nrow(huc12_target),
            paste(target_huc4s, collapse = ", ")))

# =============================================================================
# 2. BUILD BOUNDING BOX COVERING ALL TARGET HUC4s
# =============================================================================

PRISM_DIR        <- "PRISM/prism_archive"
TARGET_CACHE_DIR <- "PRISM/prism_conus"   # cropped regional tiles go here

dir.create(TARGET_CACHE_DIR, showWarnings = FALSE)

target_bbox <- huc12_target %>%
  st_bbox() %>%
  { ext(.$xmin - 0.1, .$xmax + 0.1, .$ymin - 0.1, .$ymax + 0.1) }

cat(sprintf("Bounding box (buffered):\n  xmin=%.3f  xmax=%.3f  ymin=%.3f  ymax=%.3f\n",
            target_bbox$xmin, target_bbox$xmax,
            target_bbox$ymin, target_bbox$ymax))

# =============================================================================
# 3. INDEX EXISTING PRISM ARCHIVE
# =============================================================================

prism_set_dl_dir(PRISM_DIR)

prism_files <- list.files(
  PRISM_DIR,
  pattern   = "\\.bil$",
  recursive = TRUE,
  full.names = TRUE
)

prism_dates <- str_extract(basename(prism_files), "\\d{8}") %>%
  as.Date(format = "%Y%m%d")

prism_index <- data.table(file = prism_files, date = prism_dates) %>%
  filter(!is.na(date), date >= as.Date("2009-01-01")) %>%
  arrange(date)
setDT(prism_index)

cat(sprintf("Found %d PRISM files: %s to %s\n",
            nrow(prism_index),
            min(prism_index$date),
            max(prism_index$date)))
## ----



# =============================================================================
# 4. CROP NATIONAL TILES TO TARGET BOUNDING BOX
# =============================================================================
# Output: prism_ca_hucs/PRISM_ppt_4kmD2_YYYYMMDD_CAHUC4s.tif
# Skips already-cropped files

already_cropped <- list.files(TARGET_CACHE_DIR, pattern = "\\.tif$")

to_process <- prism_index[
  !sprintf("prism_ppt_us_25m_%s_conus.tif",
           format(prism_index$date, "%Y%m%d")) %in% already_cropped
]

cat(sprintf("Files to crop: %d (skipping %d already done)\n",
            nrow(to_process),
            nrow(prism_index) - nrow(to_process)))

handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total files | ETA: :eta",
  width    = 80,
  complete = "="
))
# Convert to plain vector BEFORE the future_lapply
target_bbox_vec <- as.vector(target_bbox)  # c(xmin, xmax, ymin, ymax)

plan(multisession, workers = 10)  # start conservative
crop_one <- function(i, files, dates, target_bbox_vec, TARGET_CACHE_DIR, p) {
  on.exit(p())
  
  target_bbox <- terra::ext(target_bbox_vec)  # reconstruct inside worker
  
  r_in <- tryCatch(
    terra::rast(files[i]),
    error = function(e) {
      message(sprintf("\n  Skipping %s -- could not read: %s", files[i], e$message))
      return(NULL)
    }
  )
  if (is.null(r_in)) return(invisible(NULL))
  
  r_crop   <- terra::crop(r_in, target_bbox)
  date_str <- format(dates[i], "%Y%m%d")
  out_file <- file.path(TARGET_CACHE_DIR,
                        sprintf("prism_ppt_us_25m_%s_conus.tif", date_str))
  
  terra::writeRaster(r_crop, out_file, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  invisible(out_file)
}

with_progress({
  p <- progressor(steps = nrow(to_process))
  future_lapply(
    seq_len(nrow(to_process)),
    crop_one,
    files            = to_process$file,
    dates            = to_process$date,
    target_bbox_vec  = target_bbox_vec,
    TARGET_CACHE_DIR = TARGET_CACHE_DIR,
    p                = p,
    future.seed      = TRUE
  )
})
plan(sequential)
cat(sprintf("\nDone. Cropped tiles saved to %s/\n", TARGET_CACHE_DIR))

# =============================================================================
# 5. EXTRACT HUC12-LEVEL DAILY PRECIP FROM CROPPED TILES
# =============================================================================

tif_files <- list.files(TARGET_CACHE_DIR, pattern = "\\.tif$", full.names = TRUE)
dates     <- as.Date(str_extract(basename(tif_files), "\\d{8}"), "%Y%m%d")

tif_dt <- data.table(file = tif_files, date = dates)
tif_dt[, yr_int := data.table::year(date)]
years <- sort(unique(tif_dt$yr_int))

PRECIP_HUC12_DIR <- "PRISM/precip_huc12_tiles_conus"
dir.create(PRECIP_HUC12_DIR, showWarnings = FALSE, recursive = TRUE)

already_done <- list.files(PRECIP_HUC12_DIR, pattern = "\\.fst$")
# Split HUC12s into smaller chunks so futures do NOT receive the full 1.8 GiB object
# Easiest split: by HUC4
huc12_chunks <- split(
  huc12_target[, c("huc12", "geometry")],
  substr(huc12_target$huc12, 1, 4)
)

process_year_chunk <- function(yr, huc_chunk, tif_dt, out_dir, already_done) {
  yr <- as.integer(yr)
  yr_rows <- tif_dt[tif_dt$yr_int == yr, ]
  
  if (nrow(yr_rows) == 0) return(invisible(NULL))
  
  huc12_todo_idx <- which(
    !sprintf("precip_%s_%d.fst", huc_chunk$huc12, yr) %in% already_done
  )
  
  if (length(huc12_todo_idx) == 0) return(invisible(NULL))
  
  r_yr <- terra::rast(yr_rows$file)
  
  # Reproject only this chunk to match raster CRS
  huc12_proj <- sf::st_transform(huc_chunk, crs = terra::crs(r_yr))
  
  wide <- exactextractr::exact_extract(
    r_yr,
    huc12_proj[huc12_todo_idx, ],
    fun = "mean",
    progress = FALSE
  )
  
  for (j in seq_along(huc12_todo_idx)) {
    i <- huc12_todo_idx[j]
    huc_id <- huc_chunk$huc12[i]
    
    out_file <- file.path(out_dir, sprintf("precip_%s_%d.fst", huc_id, yr))
    
    result <- data.table::data.table(
      huc12     = huc_id,
      date      = yr_rows$date,
      precip_mm = as.numeric(wide[j, ])
    )
    
    fst::write_fst(result, out_file, compress = 75)
  }
  
  invisible(NULL)
}

process_year_chunk_safe <- function(yr, huc_chunk, tif_dt, out_dir, already_done) {
  tryCatch(
    process_year_chunk(yr, huc_chunk, tif_dt, out_dir, already_done),
    error = function(e) {
      list(
        yr       = as.integer(yr),
        chunk    = huc_chunk$huc12[1],
        error    = conditionMessage(e)
      )
    }
  )
}
# Build year × chunk tasks
tasks <- expand.grid(
  year = years,
  chunk_id = seq_along(huc12_chunks),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

plan(multisession, workers = 7)

handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total tasks | ETA: :eta",
  width    = 80,
  complete = "="
))
with_progress({
  p <- progressor(steps = nrow(tasks))
  
  results <- future_lapply(
    seq_len(nrow(tasks)),
    function(k, tasks, huc12_chunks, tif_dt, out_dir, already_done, p) {
      on.exit(p(), add = TRUE)
      yr    <- tasks$year[k]
      chunk <- huc12_chunks[[tasks$chunk_id[k]]]
      process_year_chunk_safe(yr, huc_chunk = chunk, tif_dt, out_dir, already_done)
    },
    tasks        = tasks,
    huc12_chunks = huc12_chunks,
    tif_dt       = tif_dt,
    out_dir      = PRECIP_HUC12_DIR,
    already_done = already_done,
    p            = p,
    future.seed  = TRUE
  )
})

skipped <- Filter(Negate(is.null), results)
if (length(skipped) > 0) {
  message(length(skipped), " tasks failed:")
  print(data.table::rbindlist(skipped))
} else {
  message("All tasks completed successfully.")
}


skipped_dt <- data.table::rbindlist(skipped)

missing_dt <- skipped_dt[grepl("file does not exist", error)]
crop_dt    <- skipped_dt[grepl("\\[crop\\] too few values for writing", error)]

#376 HUC12-year pairs missing due to missing data.

# =============================================================================
# 6. READ BACK AND COMBINE ALL FST FILES
# =============================================================================

files      <- list.files(PRECIP_HUC12_DIR, pattern = "\\.fst$", full.names = TRUE)
library(data.table)
library(fst)

files <- list.files(PRECIP_HUC12_DIR, pattern = "\\.fst$", full.names = TRUE)

batch_size <- 5000   
n_batches  <- ceiling(length(files) / batch_size)

precip_list <- vector("list", n_batches)

for (i in seq_len(n_batches)) {
  idx <- ((i - 1) * batch_size + 1):min(i * batch_size, length(files))
  batch_files <- files[idx]
  
  message("Reading batch ", i, " of ", n_batches)
  
  batch_dt <- rbindlist(lapply(batch_files, function(f) {
    tryCatch(
      as.data.table(read_fst(f)),
      error = function(e) {
        warning("Skipping: ", f, " — ", conditionMessage(e))
        NULL
      }
    )
  }), fill = TRUE)
  
  precip_list[[i]] <- batch_dt
  
  gc()
}

precip_all <- rbindlist(precip_list, fill = TRUE)
# Add huc4 identifier for convenience
precip_all[, huc4 := str_sub(huc12, 1, 4)]

saveRDS(precip_all, "data/precip_huc12s_conus.rds")









