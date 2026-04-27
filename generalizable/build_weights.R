# =============================================================
# build_weights.R
#
# Builds HUC12 x tract x event-day weights table.
# Policy weight is PREDETERMINED: mean active policies per tract
# across 2009 event days only, to avoid endogeneity from
# post-treatment land cover changes or flood-driven policy lapses.
#
# Columns in output:
#   huc4, upstream_huc12, tract, date
#   dist_km            -- straight-line distance (kept for reference)
#   decay_weight       -- exp(-lambda * dist_km)
#   area_sqkm_huc12    -- HUC12 area
#   policies_2009      -- mean active policies in tract across 2009 events
#   coverage_2009      -- mean active coverage in tract across 2009 events
#   w_adp              -- area * decay * policies_2009
#
# Output: one RDS per HUC4 in data/htdt_weights_by_huc4/
# =============================================================

library(data.table)
library(stringr)
library(bit64)
library(lubridate)

setwd("L:/Wetland Flood Mitigation/Paper_NFIP")

# -----------------------------
# User inputs
# -----------------------------
target_huc4s    <- "1802"
half_life_km    <- 25
distance_dir    <- "data/htdt_distance_by_huc4"
out_dir         <- "data/htdt_weights_by_huc4"
lc_path         <- "sac/data/lc_huc12_sac.rds"
policies_path   <- "data/NFIP/policies.csv"
events_path     <- "sac/data/eventdays_sac.rds"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lambda <- log(2) / half_life_km
message(sprintf("Distance decay lambda: %.4f (half-life = %d km)", lambda, half_life_km))

# -----------------------------
# Load event dates
# -----------------------------
message("Loading event dates...")
event_days <- readRDS(events_path) |> as.data.table()
event_days[, date := as.Date(date)]

event_dates_all  <- unique(event_days$date)
event_dates_2009 <- unique(event_days$date[year(event_days$date) == 2009])

cat(sprintf("All event window days: %d (%s to %s)\n",
            length(event_dates_all),
            min(event_dates_all), max(event_dates_all)))
cat(sprintf("2009 event days only:  %d\n", length(event_dates_2009)))

# -----------------------------
# Load HUC12 areas
# -----------------------------
message("Loading HUC12 area from land cover...")
lc <- readRDS(lc_path) |> as.data.table()
lc[, huc12 := str_pad(as.character(huc12), width = 12, pad = "0")]
huc12_area <- unique(lc[, .(huc12, area_sqkm_huc12)])
cat(sprintf("HUC12s with area data: %d\n", nrow(huc12_area)))

# -----------------------------
# Load policies and compute predetermined 2009 weight
# -----------------------------
message("Loading raw NFIP policies...")
policies_raw <- fread(policies_path, select = c(
  "censusTract",
  "policyEffectiveDate",
  "policyTerminationDate",
  "totalBuildingInsuranceCoverage",
  "totalContentsInsuranceCoverage"
))

policies_raw[, tract := str_pad(as.character(censusTract), width = 11, pad = "0")]
policies_raw[, `:=`(
  policyEffectiveDate   = as.Date(policyEffectiveDate),
  policyTerminationDate = as.Date(policyTerminationDate),
  total_coverage        = totalBuildingInsuranceCoverage + totalContentsInsuranceCoverage
)]
policies_raw[, c("censusTract",
                 "totalBuildingInsuranceCoverage",
                 "totalContentsInsuranceCoverage") := NULL]

# Filter to policies that could overlap 2009 event days
policies_2009_raw <- policies_raw[
  policyEffectiveDate   <= max(event_dates_2009) &
    policyTerminationDate >= min(event_dates_2009)
]
cat(sprintf("Policies overlapping 2009 event dates: %d\n", nrow(policies_2009_raw)))

# Expand to tract x 2009 event day
message("Expanding policies to 2009 event days...")
event_dt_2009 <- data.table(date = event_dates_2009)

tract_day_2009 <- policies_2009_raw[
  event_dt_2009,
  on      = .(policyEffectiveDate <= date, policyTerminationDate >= date),
  .(tract, date = i.date, total_coverage),
  nomatch = 0
][
  , .(
    n_active_policies = .N,
    total_coverage    = sum(total_coverage, na.rm = TRUE)
  ),
  by = .(tract, date)
]

ref_date <- as.Date("2009-01-01")

policies_2009_weight <- policies_raw[
  policyEffectiveDate   <= as.Date("2009-12-31") &
    policyTerminationDate >= as.Date("2009-01-01"),
  .(
    policies_2009 = .N,
    coverage_2009 = sum(total_coverage, na.rm = TRUE)
  ),
  by = tract
]

cat(sprintf("Tracts with 2009 policy weight: %d\n", nrow(policies_2009_weight)))
cat(sprintf("policies_2009 range: [%.1f, %.1f]\n",
            min(policies_2009_weight$policies_2009),
            max(policies_2009_weight$policies_2009)))

saveRDS(policies_2009_weight, "sac/data/policies_2009_weight.rds")

# -----------------------------
# Build event_id map
# -----------------------------
event_id_map <- unique(event_days[, .(date, event_id)])

# -----------------------------
# Build per HUC4
# -----------------------------
results <- lapply(target_huc4s, function(huc4_id) {
  
  message("\n=== Processing HUC4: ", huc4_id, " ===")
  
  dist_file <- file.path(distance_dir, paste0("htdt_dist_", huc4_id, ".csv"))
  if (!file.exists(dist_file)) stop("Distance file not found: ", dist_file)
  
  htdt_dist <- fread(dist_file)
  htdt_dist[, upstream_huc12 := str_pad(as.character(as.integer64(upstream_huc12)), 12, pad = "0")]
  htdt_dist[, tract := str_pad(as.character(tract), 11, pad = "0")]
  htdt_dist[, huc4  := as.character(huc4)]
  htdt_dist <- htdt_dist[huc4 == huc4_id]
  
  cat(sprintf("  HUC12-tract pairs: %d\n",        nrow(htdt_dist)))
  cat(sprintf("  Unique upstream HUC12s: %d\n",   uniqueN(htdt_dist$upstream_huc12)))
  cat(sprintf("  Unique downstream tracts: %d\n", uniqueN(htdt_dist$tract)))
  
  # ------------------------------------------------------------------
  # Static pair table: distance decay + area
  # ------------------------------------------------------------------
  weights_static <- htdt_dist[, .(upstream_huc12, tract, huc4, dist_km)]
  weights_static[, decay_weight := exp(-lambda * dist_km)]
  
  weights_static <- merge(
    weights_static,
    huc12_area,
    by.x  = "upstream_huc12",
    by.y  = "huc12",
    all.x = TRUE
  )
  
  n_missing_area <- sum(is.na(weights_static$area_sqkm_huc12))
  if (n_missing_area > 0) {
    warning(sprintf("  %d pairs missing HUC12 area -- filling with median", n_missing_area))
    weights_static[is.na(area_sqkm_huc12),
                   area_sqkm_huc12 := median(weights_static$area_sqkm_huc12, na.rm = TRUE)]
  }
  
  # ------------------------------------------------------------------
  # Join predetermined 2009 policy weight (static, not date-varying)
  # ------------------------------------------------------------------
  weights_static <- merge(
    weights_static,
    policies_2009_weight[, .(tract, policies_2009, coverage_2009)],
    by    = "tract",
    all.x = TRUE
  )
  
  weights_static[is.na(policies_2009), `:=`(policies_2009 = 0, coverage_2009 = 0)]
  
  # Combined static weight
  weights_static[, w_adp := area_sqkm_huc12 * decay_weight * policies_2009]
  
  # ------------------------------------------------------------------
  # Expand to all event dates
  # ------------------------------------------------------------------
  message("  Expanding to all event dates...")
  event_dt_all <- data.table(date = event_dates_all, dummy = 1L)
  weights_static[, dummy := 1L]
  
  weights <- merge(weights_static, event_dt_all, by = "dummy", allow.cartesian = TRUE)
  weights[, dummy := NULL]
  weights_static[, dummy := NULL]
  
  # Join event_id
  weights <- merge(weights, event_id_map, by = "date", all.x = TRUE)
  
  # ------------------------------------------------------------------
  # Diagnostics
  # ------------------------------------------------------------------
  cat(sprintf("  HUC12 x tract x day rows: %d\n",   nrow(weights)))
  cat(sprintf("  Unique event days: %d\n",           uniqueN(weights$date)))
  cat(sprintf("  dist_km range:           [%.1f, %.1f]\n",
              min(weights$dist_km),          max(weights$dist_km)))
  cat(sprintf("  decay_weight range:      [%.4f, %.4f]\n",
              min(weights$decay_weight),     max(weights$decay_weight)))
  cat(sprintf("  area_sqkm_huc12 range:   [%.2f, %.2f]\n",
              min(weights$area_sqkm_huc12),  max(weights$area_sqkm_huc12)))
  cat(sprintf("  policies_2009 range:     [%.1f, %.1f]\n",
              min(weights$policies_2009),    max(weights$policies_2009)))
  cat(sprintf("  w_adp range:             [%.4f, %.4f]\n",
              min(weights$w_adp),            max(weights$w_adp)))
  
  # ------------------------------------------------------------------
  # Save
  # ------------------------------------------------------------------
  setcolorder(weights, c(
    "huc4", "upstream_huc12", "tract", "date", "event_id",
    "dist_km", "decay_weight", "area_sqkm_huc12",
    "policies_2009", "coverage_2009", "w_adp"
  ))
  setorder(weights, upstream_huc12, tract, date)
  
  out_file <- file.path(out_dir, paste0("htdt_weights_", huc4_id, ".rds"))
  saveRDS(weights, out_file)
  message("  Saved: ", out_file)
  
  invisible(weights)
})

message("\nDone. Weight tables saved to: ", out_dir)
