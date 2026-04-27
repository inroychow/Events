
# =============================================================

library(data.table)
library(stringr)

setwd("L:/Wetland Flood Mitigation/Paper_NFIP")

# -----------------------------
# User inputs
# -----------------------------
target_huc4s <- "1802"
weights_dir  <- "sac/data"
lc_path      <- "sac/data/lc_huc12_sac.rds"
runoff_path  <- "sac/data/runoff_huc12_sac.rds"
out_dir      <- "sac/data"

# NLCD snapshot years
lc_breaks <- c(2009L, 2012L, 2015L, 2018L, 2021L, 2024L)

assign_lc_year <- function(dates) {
  yr <- as.integer(format(dates, "%Y"))
  vapply(yr, function(y) {
    valid <- lc_breaks[lc_breaks <= y]
    if (length(valid) == 0L) lc_breaks[1L] else max(valid)
  }, integer(1L))
}

lc_share_cols <- c(
  "agriculture_share", "barren_share", "developed_share",
  "forest_share", "shrub_grass_share", "water_share", "wetland_share",
  "undeveloped_share"
)

runoff_cols <- c("runoff_in", "precip_mm", "roll3_huc12", "precip_3day_in")

# =============================================================
# Load and reshape land cover (wide -> long)
# =============================================================
message("Loading and reshaping land cover shares...")
lc <- readRDS(lc_path) |> as.data.table()
lc[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]

# identify share columns -- pattern: <cover>_share_<year>
share_cols <- grep("_share_\\d{4}$", names(lc), value = TRUE)

# melt to long
lc_long <- melt(
  lc[, c("huc12", share_cols), with = FALSE],
  id.vars      = "huc12",
  measure.vars = share_cols,
  variable.name = "col_name",
  value.name    = "share_value"
)

# parse cover type and year
lc_long[, lc_year    := as.integer(str_extract(col_name, "\\d{4}$"))]
lc_long[, cover_type := str_remove(col_name, "_\\d{4}$")]
lc_long[, col_name   := NULL]

# cast back to one column per cover type, one row per huc12 x lc_year
lc <- dcast(lc_long, huc12 + lc_year ~ cover_type, value.var = "share_value")

cat(sprintf("LC rows: %d | HUC12s: %d | Years: %s\n",
            nrow(lc),
            uniqueN(lc$huc12),
            paste(sort(unique(lc$lc_year)), collapse = ", ")))

# =============================================================
# Load runoff
# =============================================================
message("Loading runoff...")
runoff <- readRDS(runoff_path) |> as.data.table()
runoff[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
runoff[, date  := as.Date(date)]

runoff_cols_keep <- c("huc12", "date", intersect(runoff_cols, names(runoff)))
runoff <- unique(runoff[, ..runoff_cols_keep])

target_huc4s = "1802"
# =============================================================
# Loop over HUC4s
# =============================================================
for (huc4_id in target_huc4s) {
  message("\n=== HUC4: ", huc4_id, " ===")
  
  # ------------------------------------------------------------------
  # Load weights (huc12 x tract x date, with all components)
  # ------------------------------------------------------------------
  w <- readRDS(file.path(weights_dir,
                         paste0("htdt_weights_", huc4_id, ".rds")))
  w[, upstream_huc12 := str_pad(as.character(upstream_huc12), 12, pad = "0")]
  w[, tract          := str_pad(as.character(tract),          11, pad = "0")]
  w[, date           := as.Date(date)]
  
  cat(sprintf("  Weight rows (huc12 x tract x date): %d\n",  nrow(w)))
  cat(sprintf("  Unique HUC12s: %d | Tracts: %d | Dates: %d\n",
              uniqueN(w$upstream_huc12),
              uniqueN(w$tract),
              uniqueN(w$date)))
  
  # ── EXPAND TO FULL EVENT WINDOWS ──────────────────────────────────
  T_pre  <- 14L
  T_post <- 14L
  
  event_dates <- unique(w$date)
  
  window_dates <- rbindlist(lapply(event_dates, function(ed) {
    data.table(event_date = ed,
               date       = seq(ed - T_pre, ed + T_post, by = "day"))
  }))
  
  w <- merge(w, window_dates,
             by.x = "date", by.y = "event_date",
             allow.cartesian = TRUE)
  w[, date := date.y][, date.y := NULL]
  w <- unique(w)
  # ──────────────────────────────────────────────────────────────────
  
  # assign lc_year to each event date
  date_lut <- unique(w[, .(date)])
  date_lut[, year_tmp := year(date)]
  date_lut[, lc_year := lc_breaks[pmax(1L, findInterval(year_tmp, lc_breaks))]]
  date_lut[, year_tmp := NULL]
  
  w <- date_lut[w, on = "date"]  
  w[, w_adp := area_sqkm_huc12 * decay_weight * policies_2009]# adp weight: area x decay x policies active downstream on this date
  
  # ================================================================
  # PART 1 — Upstream land cover (tract x date)
  # Scheme: area x decay x n_active_policies  (adp)
  # ================================================================
  message("  Building upstream land cover (area x decay x policies)...")
  
  lc_huc4 <- lc[huc12 %in% unique(w$upstream_huc12)]
  
  # join LC shares onto weight pairs by huc12 x lc_year
  lc_weighted <- merge(
    w[, .(upstream_huc12, tract, huc4, date, lc_year, w_adp)],
    lc_huc4,
    by.x  = c("upstream_huc12", "lc_year"),
    by.y  = c("huc12",          "lc_year"),
    all.x = TRUE
  )
  
  n_missing_lc <- lc_weighted[is.na(agriculture_share), .N]
  if (n_missing_lc > 0L)
    warning(sprintf("  %d huc12-lc_year pairs missing LC shares", n_missing_lc))
  
  lc_share_present <- intersect(lc_share_cols, names(lc_weighted))
  
  # normalised weighted mean per tract x date
  lc_tract <- lc_weighted[
    w_adp > 0,
    {
      denom <- sum(w_adp, na.rm = TRUE)
      out   <- lapply(lc_share_present, function(v) {
        sum(.SD[[v]] * w_adp, na.rm = TRUE) / denom
      })
      names(out) <- lc_share_present
      c(list(sum_w_adp = denom), out)
    },
    by      = .(tract, huc4, date, lc_year),
    .SDcols = lc_share_present
  ]
  
  setorder(lc_tract, tract, date)
  
  cat(sprintf("  LC output rows: %d | Tracts: %d | Dates: %d\n",
              nrow(lc_tract),
              uniqueN(lc_tract$tract),
              uniqueN(lc_tract$date)))
  cat(sprintf("  lc_year values: %s\n",
              paste(sort(unique(lc_tract$lc_year)), collapse = ", ")))
  # 
  # out_lc <- file.path(out_dir,
  #                     paste0("tract_day_lc_upstream_adp_", huc4_id, ".rds"))
  # saveRDS(lc_tract, out_lc)
  # message("  Saved LC → ", out_lc)
  
  # ================================================================
  # PART 2 — Upstream runoff (tract x date)

  # ================================================================
  message("  Building upstream runoff (policy-weighted)...")
  
  runoff_huc4 <- runoff[huc12 %in% unique(w$upstream_huc12)]
  
  # join runoff onto weight pairs by huc12 x date
  runoff_weighted <- merge(
    w[, .(upstream_huc12, tract, huc4, date, w_adp)],
    runoff_huc4,
    by.x = c("upstream_huc12", "date"),
    by.y = c("huc12",          "date"),
    all.x = TRUE
  )
  # fill missing runoff with 0
  runoff_cols_present <- intersect(runoff_cols, names(runoff_weighted))
  for (col in runoff_cols_present) {
    runoff_weighted[is.na(get(col)), (col) := 0]
  }
  
  # normalised policy-weighted mean per tract x date
  runoff_tract <- runoff_weighted[
    w_adp > 0,
    {
      denom <- sum(w_adp, na.rm = TRUE)
      out <- lapply(runoff_cols_present, function(v) {
        sum(.SD[[v]] * w_adp, na.rm = TRUE) / denom
      })
      names(out) <- paste0("upstream_", runoff_cols_present)
      out
    },
    by      = .(tract, huc4, date),
    .SDcols = runoff_cols_present
  ]
  
  setorder(runoff_tract, tract, date)
  
  cat(sprintf("  Runoff output rows: %d | Tracts: %d | Dates: %d\n",
              nrow(runoff_tract),
              uniqueN(runoff_tract$tract),
              uniqueN(runoff_tract$date)))
  cat(sprintf("  Date range: %s to %s\n",
              min(runoff_tract$date),
              max(runoff_tract$date)))
  tract_day <- merge(
    lc_tract,
    runoff_tract,
    by  = c("tract", "huc4", "date"),
    all = TRUE
  )
  
  out_path <- file.path(out_dir, paste0("weighted_lc_runoff_", huc4_id, ".rds"))
  saveRDS(tract_day, out_path)
  message("  Saved combined → ", out_path)
}

message("\nDone.")
