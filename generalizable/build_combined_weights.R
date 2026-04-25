# =============================================================
# build_combined_weights.R
#
# Builds HUC12-tract combined weights incorporating:
#   1. Distance decay  -- exp(-lambda * dist_km)
#   2. HUC12 area      -- area_sqkm_huc12
#   3. 2009 tract policies -- predetermined insured exposure
#
# combined_weight = exp(-lambda * dist_km) * area_sqkm_huc12 * policies_2009
#
# Output: htdt_dist files, one per HUC4, with combined_weight column appended
#         saved to data/htdt_distance_by_huc4_weighted/
#
# These weighted distance files replace the unweighted ones used in
# build_tract_day_runoff_expdecay() and build_tract_year_lc_expdecay()
# in script 2. Downstream aggregations then use combined_weight in place
# of decay_weight, so no additional weighting is needed in script 4.
# =============================================================

library(data.table)
library(stringr)
library(bit64)
library(dplyr)

setwd("L:/Wetland Flood Mitigation/Paper_NFIP")

# -----------------------------
# User inputs
# -----------------------------
target_huc4s   <- "1802"           # can be a vector e.g. c("1802", "1801")
half_life_km   <- 25               # distance decay half-life in km
distance_dir   <- "data/htdt_distance_by_huc4"
out_dir        <- "data/htdt_distance_by_huc4_weighted"
panel2_path    <- "sac/data/panel2_sac.rds"
lc_path        <- "sac/data/lc_huc12_sac.rds"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load inputs
# -----------------------------
message("Loading panel2 for 2009 policy weights...")
panel2 <- readRDS(panel2_path) |> as.data.table()
panel2[, tract := str_pad(as.character(tract), width = 11, side = "left", pad = "0")]

message("Loading HUC12 land cover for area...")
lc <- readRDS(lc_path) |> as.data.table()
lc[, huc12 := str_pad(as.character(huc12), width = 12, side = "left", pad = "0")]

# -----------------------------
# Build 2009 tract policy weights
# -----------------------------
message("Building 2009 tract policy weights...")
tract_policies_2009 <- panel2[
  year == 2009,
  .(policies_2009 = sum(active_policies, na.rm = TRUE)),
  by = tract
]

cat(sprintf("  Tracts with 2009 policy data: %d\n", nrow(tract_policies_2009)))
cat(sprintf("  Tracts with zero 2009 policies: %d\n",
            sum(tract_policies_2009$policies_2009 == 0)))

# -----------------------------
# Build HUC12 area lookup
# -----------------------------
huc12_area <- unique(lc[, .(huc12, area_sqkm_huc12)])
cat(sprintf("  HUC12s with area data: %d\n", nrow(huc12_area)))

# -----------------------------
# Build combined weights per HUC4
# -----------------------------
lambda <- log(2) / half_life_km
message(sprintf("Distance decay lambda: %.4f (half-life = %d km)", lambda, half_life_km))

results <- lapply(target_huc4s, function(huc4_id) {
  
  message("\nProcessing HUC4: ", huc4_id)
  
  # load distance file
  dist_file <- file.path(distance_dir, paste0("htdt_dist_", huc4_id, ".csv"))
  if (!file.exists(dist_file)) stop("Distance file not found: ", dist_file)
  
  htdt_dist <- fread(dist_file)
  
  # standardize IDs
  htdt_dist[, upstream_huc12 := str_pad(
    as.character(as.integer64(upstream_huc12)), 12, pad = "0"
  )]
  htdt_dist[, tract := str_pad(as.character(tract), 11, pad = "0")]
  htdt_dist[, huc4  := as.character(huc4)]
  htdt_dist <- htdt_dist[huc4 == huc4_id]
  
  cat(sprintf("  HUC12-tract pairs: %d\n", nrow(htdt_dist)))
  cat(sprintf("  Unique upstream HUC12s: %d\n", uniqueN(htdt_dist$upstream_huc12)))
  cat(sprintf("  Unique downstream tracts: %d\n", uniqueN(htdt_dist$tract)))
  
  # -----------------------------
  # Component 1: distance decay
  # -----------------------------
  htdt_dist[, decay_weight := exp(-lambda * dist_km)]
  
  # -----------------------------
  # Component 2: HUC12 area
  # -----------------------------
  htdt_dist <- merge(
    htdt_dist,
    huc12_area,
    by.x = "upstream_huc12",
    by.y = "huc12",
    all.x = TRUE
  )
  
  n_missing_area <- sum(is.na(htdt_dist$area_sqkm_huc12))
  if (n_missing_area > 0) {
    warning(sprintf("  %d HUC12-tract pairs missing area — filling with median area", n_missing_area))
    median_area <- median(htdt_dist$area_sqkm_huc12, na.rm = TRUE)
    htdt_dist[is.na(area_sqkm_huc12), area_sqkm_huc12 := median_area]
  }
  
  # -----------------------------
  # Component 3: 2009 tract policy weight
  # -----------------------------
  htdt_dist <- merge(
    htdt_dist,
    tract_policies_2009,
    by = "tract",
    all.x = TRUE
  )
  
  n_missing_pol <- sum(is.na(htdt_dist$policies_2009))
  if (n_missing_pol > 0) {
    warning(sprintf("  %d tracts missing 2009 policies — filling with 0", n_missing_pol))
    htdt_dist[is.na(policies_2009), policies_2009 := 0]
  }
  
  # -----------------------------
  # Combined weight
  # -----------------------------
  htdt_dist[, combined_weight := decay_weight * area_sqkm_huc12 * policies_2009]
  
  # summary diagnostics
  cat(sprintf("  Combined weight range: [%.4f, %.4f]\n",
              min(htdt_dist$combined_weight),
              max(htdt_dist$combined_weight)))
  cat(sprintf("  HUC12-tract pairs with zero combined weight: %d (%.1f%%)\n",
              sum(htdt_dist$combined_weight == 0),
              100 * mean(htdt_dist$combined_weight == 0)))
  
  # -----------------------------
  # Save
  # -----------------------------
  out_file <- file.path(out_dir, paste0("htdt_dist_weighted_", huc4_id, ".csv"))
  fwrite(htdt_dist, out_file)
  message("  Saved: ", out_file)
  
  invisible(htdt_dist)
})

message("\nDone. Weighted distance files saved to: ", out_dir)
message("Update distance_dir in build_tract_day_runoff_expdecay() and")
message("build_tract_year_lc_expdecay() to point to: ", out_dir)
message("And replace w = decay_weight with w = combined_weight in both functions.")
