library(sf)
library(data.table)
library(dplyr)
library(stringr)
library(bit64)
library(tidyr)

# ------------------------------------------------------------
# Inputs needed:
# precip_event_huc12: huc12-date runoff data
# htdt_distance_by_huc4_weighted/: weighted distance files from build_combined_weights.R
# lc_huc12_sac.rds: HUC12 land cover shares
# ------------------------------------------------------------

target_huc4s <- "1802"

# =============================================================
# RUNOFF
# build_tract_day_runoff_expdecay()
#
# Aggregates HUC12-day runoff to tract-day using combined weights
# (distance decay x HUC12 area x 2009 tract policies) from
# build_combined_weights.R. No weighting recomputed here.
# =============================================================

build_tract_day_runoff_expdecay <- function(
    huc4_id,
    precip_event_huc12,
    distance_dir = "data/htdt_distance_by_huc4_weighted",
    out_dir      = "sac/data"
) {
  
  message("Processing HUC4: ", huc4_id)
  
  # -----------------------------
  # Load combined-weight distance file
  # -----------------------------
  dist_file <- file.path(distance_dir, paste0("htdt_dist_weighted_", huc4_id, ".csv"))
  
  if (!file.exists(dist_file)) {
    stop("Weighted distance file not found: ", dist_file,
         "\nRun build_combined_weights.R first.")
  }
  
  htdt_dist <- fread(dist_file)
  
  htdt_dist[, upstream_huc12 := str_pad(
    as.character(as.integer64(upstream_huc12)), 12, pad = "0"
  )]
  htdt_dist[, tract := str_pad(as.character(tract), 11, pad = "0")]
  htdt_dist[, huc4  := as.character(huc4)]
  htdt_dist <- htdt_dist[huc4 == huc4_id]
  
  # -----------------------------
  # Prep runoff data
  # -----------------------------
  runoff_dt <- as.data.table(precip_event_huc12)
  
  runoff_dt[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
  runoff_dt[, huc4  := str_sub(huc12, 1, 4)]
  runoff_dt[, date  := as.Date(date)]
  
  runoff_dt <- runoff_dt[huc4 == huc4_id]
  
  # -----------------------------
  # Join HUC12 runoff to downstream tracts
  # -----------------------------
  tract_runoff_long <- merge(
    htdt_dist,
    runoff_dt,
    by.x = "upstream_huc12",
    by.y = "huc12",
    allow.cartesian = TRUE
  )
  
  # -----------------------------
  # Aggregate to tract-day using combined_weight
  # (distance decay x HUC12 area x 2009 policies — built in build_combined_weights.R)
  # -----------------------------
  tract_day_runoff <- tract_runoff_long[
    ,
    .(
      upstream_runoff_in      = weighted.mean(runoff_in,      w = combined_weight, na.rm = TRUE),
      upstream_precip_3day_in = weighted.mean(precip_3day_in, w = combined_weight, na.rm = TRUE)
    ),
    by = .(tract, date)
  ]
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(out_dir, paste0("tract_day_dd_runoff_", huc4_id, ".rds"))
  saveRDS(tract_day_runoff, out_file)
  message("Saved: ", out_file)
  
  return(tract_day_runoff)
}

# -----------------------------
# Run
# -----------------------------
runoff_list <- lapply(target_huc4s, function(huc4) {
  build_tract_day_runoff_expdecay(
    huc4_id           = huc4,
    precip_event_huc12 = precip_event_huc12
  )
})

# tract_day_runoff_all <- rbindlist(runoff_list, fill = TRUE)
# saveRDS(
#   tract_day_runoff_all,
#   "data/tract_day_upstream_runoff_expdecay_all_huc4s.rds"
# )


# =============================================================
# LAND COVER
# build_tract_year_lc_expdecay()
#
# Aggregates HUC12-year land cover to tract-year using combined weights
# (distance decay x HUC12 area x 2009 tract policies) from
# build_combined_weights.R. No weighting recomputed here.
# =============================================================

lc_shares <- readRDS("sac/data/lc_huc12_sac.rds")

lc_long <- lc_shares %>%
  pivot_longer(
    cols = matches("_(share|sqkm)_\\d{4}$"),
    names_to  = c("group", "metric", "year"),
    names_pattern = "(.+)_(share|sqkm)_(\\d{4})",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols    = c(huc12, huc4, year),
    names_from = c(group, metric),
    values_from = value
  ) %>%
  mutate(
    year  = as.integer(year),
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    huc4  = str_sub(huc12, 1, 4)
  )

build_tract_year_lc_expdecay <- function(
    huc4_id,
    lc_long,
    distance_dir = "data/htdt_distance_by_huc4_weighted",
    out_dir      = "sac/data"
) {
  
  message("Processing HUC4: ", huc4_id)
  
  # -----------------------------
  # Load combined-weight distance file
  # -----------------------------
  dist_file <- file.path(distance_dir, paste0("htdt_dist_weighted_", huc4_id, ".csv"))
  
  if (!file.exists(dist_file)) {
    stop("Weighted distance file not found: ", dist_file,
         "\nRun build_combined_weights.R first.")
  }
  
  htdt_dist <- fread(dist_file)
  
  htdt_dist[, upstream_huc12 := str_pad(
    as.character(as.integer64(upstream_huc12)), 12, pad = "0"
  )]
  htdt_dist[, tract := str_pad(as.character(tract), 11, pad = "0")]
  htdt_dist[, huc4  := as.character(huc4)]
  htdt_dist <- htdt_dist[huc4 == huc4_id]
  
  # -----------------------------
  # Prep land cover data
  # -----------------------------
  lc_dt <- as.data.table(lc_long)
  lc_dt[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
  lc_dt[, huc4  := str_sub(huc12, 1, 4)]
  lc_dt <- lc_dt[huc4 == huc4_id]
  
  # -----------------------------
  # Join HUC12 land cover to downstream tracts
  # -----------------------------
  tract_huc_lc <- merge(
    htdt_dist,
    lc_dt,
    by.x = "upstream_huc12",
    by.y = "huc12",
    allow.cartesian = TRUE
  )
  
  lc_vars <- c(
    "agriculture_share",
    "barren_share",
    "developed_share",
    "forest_share",
    "shrub_grass_share",
    "water_share",
    "wetland_share",
    "undeveloped_share"
  )
  
  lc_vars <- intersect(lc_vars, names(tract_huc_lc))
  
  # -----------------------------
  # Aggregate to tract-year using combined_weight
  # (distance decay x HUC12 area x 2009 policies — built in build_combined_weights.R)
  # -----------------------------
  tract_year_lc_expdecay <- tract_huc_lc[
    ,
    c(
      lapply(.SD, function(x) weighted.mean(x, w = combined_weight, na.rm = TRUE)),
      list(
        mean_upstream_dist_km = weighted.mean(dist_km,        w = combined_weight, na.rm = TRUE),
        n_upstream_huc12      = uniqueN(upstream_huc12)
      )
    ),
    by      = .(tract, year),
    .SDcols = lc_vars
  ]
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(out_dir, paste0("tract_year_lc_dd_", huc4_id, ".rds"))
  saveRDS(tract_year_lc_expdecay, out_file)
  message("Saved: ", out_file)
  
  return(tract_year_lc_expdecay)
}

# -----------------------------
# Run
# -----------------------------
tract_year_lc_expdecay_1802 <- build_tract_year_lc_expdecay(
  huc4_id  = "1802",
  lc_long  = lc_long
)
