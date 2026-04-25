setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")

library(data.table)
library(dplyr)
library(fixest)
library(lubridate)
library(tidyverse)
library(stringr)

# target_huc4s <- c("1806", "1801", "1805", "1804", "1803")

target_huc4s = "1802"
# Run build_lc, build_nfip, and build_events first
data_path <- "sac/data"

claims_all <- readRDS(file.path(data_path, "claims_sac.rds"))
policies_all <- readRDS(file.path(data_path, "policies_sac.rds"))
event_days_all <- readRDS(file.path(data_path, "eventdays_post2009_sac.rds"))
tract_year_lc <- readRDS(file.path(data_path, "tract_year_lc_sac.rds"))
# or: tract_year_lc_cahucs_aw.rds

setDT(claims_all)
setDT(policies_all)
setDT(event_days_all)
tract_year_lc_dt <- as.data.table(tract_year_lc)


# -----------------------------#
# 1. Standardize inputs
# -----------------------------#

claims_all[, tract := as.character(tract)]
policies_all[, tract := as.character(tract)]
tract_year_lc_dt[, tract := as.character(tract)]

claims_all[, dateOfLoss := as.Date(dateOfLoss)]
policies_all[, policyEffectiveDate := as.Date(policyEffectiveDate)]
policies_all[, policyTerminationDate := as.Date(policyTerminationDate)]

event_days_all[, event_start := as.Date(event_start)]
event_days_all[, event_end   := as.Date(event_end)]
event_days_all[, peak_date   := as.Date(peak_date)]

tract_year_lc_dt[, year := as.integer(year)]

library(data.table)
library(dplyr)
library(fixest)
library(lubridate)
library(tidyverse)

# =============================================================
# build_panel.R  –  generalizable panel builder
# 
# =============================================================

build_panel <- function(
    target_huc4s,           # character vector of HUC4 codes to include
    region_tag,             # short string used in file paths, e.g. "sac", "delta"
    data_dir       = NULL,  # base data dir; defaults to "{region_tag}/data"
    out_dir        = NULL,  # output dir;   defaults to data_dir
    panel_start    = as.Date("2009-01-01"),
    T_pre          = 7L,    # days before event_start for event-study window
    T_post         = 7L,    # days after  event_start for event-study window
    lc_area_weighted = FALSE, # TRUE → use tract_year_lc_aw.rds
    save_main_panel  = TRUE,
    save_event_panel = TRUE,
    runoff_file = "tract_day_dd_runoff_1802.rds"
) {
  
  # ------------------------------------------------------------------
  #   Resolve file paths
  # ------------------------------------------------------------------
  if (is.null(data_dir)) data_dir <- file.path(region_tag, "data")
  if (is.null(out_dir))  out_dir  <- data_dir
  
  # lc_file <- if (lc_area_weighted) "tract_year_lc_sac_aw.rds" else "tract_year_lc_sac.rds" # this is not dd weighted
  lc_file <- paste0("tract_year_lc_dd_", target_huc4s, ".rds") #dd weighted
  
  message("=== build_panel | region: ", region_tag,
          " | HUC4s: ", paste(target_huc4s, collapse = ", "), " ===")
  
  # ------------------------------------------------------------------
  #  Load + standardize data inputs
  # ------------------------------------------------------------------
  claims_raw   <- readRDS(file.path(data_dir, "claims_sac.rds"))        |> as.data.table()
  policies_raw <- readRDS(file.path(data_dir, "policies_sac.rds"))      |> as.data.table()
  events_raw   <- readRDS(file.path(data_dir, "eventdays_post2009_sac.rds")) |> as.data.table()
  lc_raw       <- readRDS(file.path(data_dir, lc_file))                 |> as.data.table()
  runoff_raw <- readRDS(file.path(data_dir, runoff_file)) |> as.data.table()
  
  # character IDs
  for (dt in list(claims_raw, policies_raw, lc_raw)) {
    if ("tract" %in% names(dt)) dt[, tract := as.character(tract)]
    if ("huc4"  %in% names(dt)) dt[, huc4  := as.character(huc4)]
  }
  
  # dates 
  claims_raw[, dateOfLoss            := as.Date(dateOfLoss)]
  policies_raw[, policyEffectiveDate  := as.Date(policyEffectiveDate)]
  policies_raw[, policyTerminationDate := as.Date(policyTerminationDate)]
  events_raw[, event_start := as.Date(event_start)]
  events_raw[, event_end   := as.Date(event_end)]
  events_raw[, peak_date   := as.Date(peak_date)]
  lc_raw[, year := as.integer(year)]
  
  runoff_raw[, tract := as.character(tract)]
  runoff_raw[, date := as.Date(date)]
  
  # ------------------------------------------------------------------
  # Filter to target HUC4s
  # ------------------------------------------------------------------
  filter_huc4 <- function(dt, huc4s) {
    if ("huc4" %in% names(dt)) dt[huc4 %in% huc4s] else dt
  }
  
  claims   <- filter_huc4(claims_raw,   target_huc4s)
  policies <- filter_huc4(policies_raw, target_huc4s)
  lc_dt    <- filter_huc4(lc_raw,       target_huc4s)
  # events are region-wide (precipitation-based); keep all
  events   <- copy(events_raw)
  
  message("  claims rows: ",   nrow(claims),
          " | policy rows: ",  nrow(policies),
          " | lc rows: ",      nrow(lc_dt))
  
  # ------------------------------------------------------------------
  # Aggregate claims --> tract-day
  # ------------------------------------------------------------------
  claims_day <- claims[
    , .(
      tract_claims    = sum(tract_claims,    na.rm = TRUE),
      tract_total_paid = sum(tract_total_paid, na.rm = TRUE)
    ),
    by = .(tract, date = dateOfLoss)
  ]
  
  # ------------------------------------------------------------------
  # Build full tract - day panel
  # ------------------------------------------------------------------
  all_tracts <- unique(c(claims$tract, policies$tract, lc_dt$tract))
  
  panel_end <- max(
    max(claims$dateOfLoss,               na.rm = TRUE),
    max(policies$policyTerminationDate,  na.rm = TRUE),
    max(events$event_end,                na.rm = TRUE),
    na.rm = TRUE
  )
  
  panel <- CJ(
    tract = all_tracts,
    date  = seq(panel_start, panel_end, by = "day"),
    unique = TRUE
  )
  
  # ------------------------------------------------------------------
  #  Join claims to panel object
  # ------------------------------------------------------------------
  panel <- claims_day[panel, on = .(tract, date)]
  panel[is.na(tract_claims),     tract_claims     := 0L]
  panel[is.na(tract_total_paid), tract_total_paid := 0]
  
  # ------------------------------------------------------------------
  # Get policies per tract day
  # ------------------------------------------------------------------
  policy_daily <- policies[
    , .(tract, start = policyEffectiveDate,
        end = policyTerminationDate, total_coverage)
  ][
    panel,
    .(tract, date = i.date, total_coverage),
    on  = .(tract, start <= date, end >= date),
    allow.cartesian = TRUE,
    nomatch = 0L
  ][
    , .(
      active_policies = .N,
      active_coverage = sum(total_coverage, na.rm = TRUE)
    ),
    by = .(tract, date)
  ]
  
  panel <- policy_daily[panel, on = .(tract, date)]
  panel[is.na(active_policies), active_policies := 0L]
  panel[is.na(active_coverage), active_coverage := 0]
  
  # ------------------------------------------------------------------
  # Extreme event indicators
  # ------------------------------------------------------------------
  event_daily <- unique(
    events[, .(date, extreme_event, event_id, event_start, event_end,
               peak_date, peak_precip, precip_mm, roll3, is_peak_day)]
  )
  
  panel <- event_daily[panel, on = "date"]
  panel[is.na(extreme_event), extreme_event := 0L]
  panel[is.na(is_peak_day),   is_peak_day   := 0L]
  
  
  # ------------------------------------------------------------------
  # Join tract-day upstream runoff exposure
  # ------------------------------------------------------------------
  
  panel <- runoff_raw[panel, on = .(tract, date)]
  
  panel[is.na(upstream_runoff_in), upstream_runoff_in := 0]
  panel[is.na(upstream_precip_3day_in), upstream_precip_3day_in := 0]
  # ------------------------------------------------------------------
  #  Land-cover year assignment  (most recent NLCD snapshot)
  # ------------------------------------------------------------------
  lc_breaks <- c(2009L, 2012L, 2015L, 2018L, 2021L, 2024L)
  
  assign_lc_year <- function(d) {
    yr <- as.integer(format(d, "%Y"))
    # find largest lc_break that is <= yr
    vapply(yr, function(y) max(lc_breaks[lc_breaks <= y]), integer(1L))
  }
  
  panel[, lc_year := assign_lc_year(date)]
  
  setnames(lc_dt, "year", "lc_year")
  panel <- lc_dt[panel, on = .(tract, lc_year)]
  
  # ------------------------------------------------------------------
  # Time vars and outcomes in log
  # ------------------------------------------------------------------
  panel[, `:=`(
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    yday  = as.integer(format(date, "%j")),
    dow   = weekdays(date)
  )]
  
  panel[, `:=`(
    log_total_paid      = log1p(tract_total_paid),
    log_active_coverage = log1p(active_coverage)
  )]
  
  # ------------------------------------------------------------------
  # LC 
  # ------------------------------------------------------------------
  landcover_vars <- c(
    "agriculture_share", "barren_share", "developed_share",
    "forest_share", "shrub_grass_share", "water_share", "wetland_share"
  )
  lc_present <- intersect(landcover_vars, names(panel))  # only those that exist
  
  for (v in lc_present) {
    panel[, paste0("event_x_", v) := .SD[[1L]] * extreme_event, .SDcols = v]
  }
  
  panel[, event_x_upstream_runoff := extreme_event * upstream_runoff_in]
  
  #Make outcome - fraction of total paid / active coverage
  panel[, claims_per_coverage := fifelse(
    active_coverage > 0, tract_total_paid / active_coverage, NA_real_
  )]
  panel[!is.na(event_id), tau := as.integer(date - event_start)]
  
  # ------------------------------------------------------------------
  # Save main panel
  # ------------------------------------------------------------------
  out_path_main <- file.path(out_dir, paste0("panel_", region_tag, ".rds"))
  if (save_main_panel) {
    saveRDS(panel, out_path_main)
    message("  Saved main panel → ", out_path_main)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ------------------------------------------------------------------
  # Event-study panel starts here
  # ------------------------------------------------------------------
  events_uniq <- panel[
    !is.na(event_id) & !is.na(event_start),
    .(event_start = event_start[1L], event_end = event_end[1L]),
    by = event_id
  ]
  
  event_windows <- events_uniq[, .(
    date = seq(event_start - T_pre, event_start + T_post, by = "day")
  ), by = .(event_id, event_start, event_end)]
  
  event_windows[, tau := as.integer(date - event_start)]
  
  tracts <- unique(panel[, .(tract)])
  tracts[, join_key := 1L]
  event_windows[, join_key := 1L]
  
  event_windows_tract <- tracts[event_windows, on = .(join_key), allow.cartesian = TRUE]
  
  tracts[, join_key := NULL]
  event_windows[, join_key := NULL]
  event_windows_tract[, join_key := NULL]
  
  panel2 <- merge(
    event_windows_tract,
    panel,
    by      = c("tract", "date"),
    all.x   = TRUE,
    allow.cartesian = TRUE
  )
  panel2 <- panel2[!is.na(date)]
  
  # get rid of duplicate columns from the merge
  setnames(panel2,
           old = c("event_id.x", "event_start.x", "event_end.x", "tau.x"),
           new = c("event_id",   "event_start",   "event_end",   "tau"))
  
  drop_cols <- intersect(
    c("event_id.y", "event_start.y", "event_end.y", "tau.y"),
    names(panel2)
  )
  panel2[, (drop_cols) := NULL]
  
  panel2[, tau_factor := relevel(factor(tau), ref = "-1")]
  panel2[, in_event_window := fifelse(
    date >= event_start & date <= event_end, 1L, 0L
  )]
  panel2[, event_year := year(event_start)]
  
  out_path_event <- file.path(out_dir, paste0("panel2_", region_tag, ".rds"))
  if (save_event_panel) {
    saveRDS(panel2, out_path_event)
    message("  Saved event panel → ", out_path_event)
  }
  
  invisible(list(panel = panel, panel2 = panel2))
}


# =============================================================
#calls
# =============================================================

# can evetually do this with a vector of all the HUC4s in conus using target_huc4s
build_panel(
  target_huc4s      = "1802",
  region_tag        = "sac",
  lc_area_weighted  = FALSE,
  runoff_file       = "tract_day_dd_runoff_1802.rds",
  T_pre             = 7L,
  T_post            = 7L,
  save_main_panel   = FALSE,
  save_event_panel  = TRUE
)
