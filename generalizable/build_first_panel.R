# =============================================================================
# build_panel.R
# =============================================================================

setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")

library(data.table)
library(dplyr)
library(fixest)
library(lubridate)
library(tidyverse)
library(stringr)

target_huc4s <- "1802"

build_panel <- function(
    target_huc4s,
    region_tag,
    data_dir         = "sac/data",
    out_dir          = NULL,
    panel_start      = as.Date("2009-01-01"),
    T_pre            = 7L,
    T_post           = 7L,
    save_main_panel  = TRUE,
    save_event_panel = TRUE,
    upstream_file    = NULL
) {
  
  # ------------------------------------------------------------------
  # Resolve file paths
  # ------------------------------------------------------------------
  if (is.null(data_dir))   data_dir <- file.path(region_tag, "sac/data")
  if (is.null(out_dir))    out_dir  <- data_dir
  if (is.null(upstream_file))
    upstream_file <- paste0("weighted_lc_runoff_", target_huc4s, ".rds")
  
  message("=== build_panel | region: ", region_tag,
          " | HUC4s: ", paste(target_huc4s, collapse = ", "), " ===")
  
  # ------------------------------------------------------------------
  # Load + standardize data inputs
  # ------------------------------------------------------------------
  claims_raw          <- readRDS(file.path(data_dir, "claims_sac.rds"))             |> as.data.table()
  policies_raw        <- readRDS(file.path(data_dir, "policies_sac.rds"))           |> as.data.table()
  events_raw          <- readRDS("sac/data/eventdays_sac.rds")                        |> as.data.table()
  upstream_raw        <- readRDS(file.path(data_dir, upstream_file))                |> as.data.table()
  weights_raw <- readRDS(file.path("sac/data", 
                                   paste0("htdt_weights_", target_huc4s, ".rds"))) |> 
    as.data.table()
  
  policies_2009_raw <- unique(weights_raw[, .(tract, policies_2009, coverage_2009)])
  policies_2009_raw[, tract := str_pad(as.character(tract), 11, pad = "0")]  
  # character IDs
  for (dt in list(claims_raw, policies_raw)) {
    if ("tract" %in% names(dt)) dt[, tract := as.character(tract)]
    if ("huc4"  %in% names(dt)) dt[, huc4  := as.character(huc4)]
  }
  
  upstream_raw[,      tract := str_pad(as.character(tract), 11, pad = "0")]
  upstream_raw[,      date  := as.Date(date)]
  policies_2009_raw[, tract := str_pad(as.character(tract), 11, pad = "0")]
  
  # dates
  claims_raw[,   dateOfLoss             := as.Date(dateOfLoss)]
  policies_raw[, policyEffectiveDate    := as.Date(policyEffectiveDate)]
  policies_raw[, policyTerminationDate  := as.Date(policyTerminationDate)]
  events_raw[,   event_start            := as.Date(event_start)]
  events_raw[,   event_end              := as.Date(event_end)]
  events_raw[,   peak_date              := as.Date(peak_date)]
  
  # ------------------------------------------------------------------
  # Filter to target HUC4s
  # ------------------------------------------------------------------
  filter_huc4 <- function(dt, huc4s) {
    if ("huc4" %in% names(dt)) dt[huc4 %in% huc4s] else dt
  }
  
  claims   <- filter_huc4(claims_raw,   target_huc4s)
  policies <- filter_huc4(policies_raw, target_huc4s)
  events   <- copy(events_raw)
  
  message("  claims rows: ",  nrow(claims),
          " | policy rows: ", nrow(policies))
  
  # ------------------------------------------------------------------
  # Aggregate claims --> tract-day
  # ------------------------------------------------------------------
  claims_day <- claims[
    , .(
      tract_claims     = sum(tract_claims,     na.rm = TRUE),
      tract_total_paid = sum(tract_total_paid, na.rm = TRUE)
    ),
    by = .(tract, date = dateOfLoss)
  ]
  
  # ------------------------------------------------------------------
  # Build full tract x day panel skeleton
  # ------------------------------------------------------------------
  all_tracts <- unique(c(claims$tract, policies$tract, upstream_raw$tract))
  
  panel_end <- max(
    max(claims$dateOfLoss,              na.rm = TRUE),
    max(policies$policyTerminationDate, na.rm = TRUE),
    max(events$event_end,               na.rm = TRUE),
    na.rm = TRUE
  )
  
  panel <- CJ(
    tract  = all_tracts,
    date   = seq(panel_start, panel_end, by = "day"),
    unique = TRUE
  )
  
  # ------------------------------------------------------------------
  # Join claims
  # ------------------------------------------------------------------
  panel <- claims_day[panel, on = .(tract, date)]
  panel[is.na(tract_claims),     tract_claims     := 0L]
  panel[is.na(tract_total_paid), tract_total_paid := 0]
  
  # ------------------------------------------------------------------
  # Join active policies per tract-day (time-varying, for controls)
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
  # Join predetermined 2009 policy weight (static, tract-level)
  # ------------------------------------------------------------------
  panel <- merge(panel, policies_2009_raw, by = "tract", all.x = TRUE)
  panel[is.na(policies_2009), `:=`(policies_2009 = 0, coverage_2009 = 0)]
  
  # ------------------------------------------------------------------
  # Join event indicators
  # ------------------------------------------------------------------
  event_daily <- unique(
    events[, .(date, extreme_event, event_id, event_start, event_end,
               peak_date, peak_precip, precip_mm, roll3, is_peak_day)]
  )
  
  panel <- event_daily[panel, on = "date"]
  panel[is.na(extreme_event), extreme_event := 0L]
  panel[is.na(is_peak_day),   is_peak_day   := 0L]
  
  # ------------------------------------------------------------------
  # Join upstream LC + runoff (single combined file)
  # ------------------------------------------------------------------
  panel <- upstream_raw[panel, on = .(tract, date)]
  
  # fill NAs for all upstream columns with 0
  upstream_cols <- grep("^upstream_|_share$", names(panel), value = TRUE)
  for (col in upstream_cols) {
    panel[is.na(get(col)), (col) := 0]
  }
  
  # ------------------------------------------------------------------
  # Time variables
  # ------------------------------------------------------------------
  panel[, `:=`(
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    yday  = as.integer(format(date, "%j")),
    dow   = weekdays(date)
  )]
  
  # ------------------------------------------------------------------
  # Outcome variables
  # ------------------------------------------------------------------
  panel[, `:=`(
    log_total_paid      = log1p(tract_total_paid),
    log_active_coverage = log1p(active_coverage)
  )]
  
  panel[, claims_per_coverage := fifelse(
    active_coverage > 0, tract_total_paid / active_coverage, NA_real_
  )]
  
  # panel[!is.na(event_id), tau := as.integer(date - event_start)]
  
  # # ------------------------------------------------------------------
  # # Interaction terms
  # # ------------------------------------------------------------------
  # # upstream runoff x predetermined exposure
  # if ("upstream_runoff_in" %in% names(panel))
  #   panel[, runoff_x_policies2009 := upstream_runoff_in * policies_2009]
  # 
  # # event x upstream runoff
  # if ("upstream_runoff_in" %in% names(panel))
  #   panel[, event_x_upstream_runoff := extreme_event * upstream_runoff_in]
  # 
  # # event x upstream runoff x predetermined exposure (triple interaction)
  # if ("upstream_runoff_in" %in% names(panel))
  #   panel[, event_x_runoff_x_policies2009 := extreme_event * upstream_runoff_in * policies_2009]
  # 
  # # ------------------------------------------------------------------
  # # Save main panel
  # # ------------------------------------------------------------------
  # if (save_main_panel) {
  #   out_path_main <- file.path(out_dir, paste0("panel_", region_tag, ".rds"))
  #   saveRDS(panel, out_path_main)
  #   message("  Saved main panel → ", out_path_main)
  # }
  
  # ------------------------------------------------------------------
  # Event-study panel
  # ------------------------------------------------------------------
  events_uniq <- panel[
    !is.na(event_id) & !is.na(peak_date),
    .(
      event_start = event_start[1L],
      event_end   = event_end[1L],
      peak_date   = peak_date[1L]
    ),
    by = event_id
  ]
  
  event_windows <- events_uniq[, .(
    date = seq(peak_date - T_pre, peak_date + T_post, by = "day")
  ), by = .(event_id, event_start, event_end, peak_date)]
  
  event_windows[, tau := as.integer(date - peak_date)]
  
  tracts <- unique(panel[, .(tract)])
  tracts[,        join_key := 1L]
  event_windows[, join_key := 1L]
  
  event_windows_tract <- tracts[event_windows, on = .(join_key), allow.cartesian = TRUE]
  tracts[,        join_key := NULL]
  event_windows[, join_key := NULL]
  event_windows_tract[, join_key := NULL]
  
  panel2 <- merge(
    event_windows_tract,
    panel,
    by              = c("tract", "date"),
    all.x           = TRUE,
    allow.cartesian = TRUE
  )
  panel2 <- panel2[!is.na(date)]
  
  setnames(panel2,
           old = c("event_id.x", "event_start.x", "event_end.x", "peak_date.x", "tau.x"),
           new = c("event_id",   "event_start",   "event_end",   "peak_date",   "tau"),
           skip_absent = TRUE)
  
  drop_cols <- intersect(
    c("event_id.y", "event_start.y", "event_end.y", "peak_date.y", "tau.y"),
    names(panel2)
  )
  panel2[, (drop_cols) := NULL]
  
  panel2[, tau_factor := relevel(factor(tau), ref = "-1")]
  panel2[, in_event_window := fifelse(date >= event_start & date <= event_end, 1L, 0L)]
  panel2[, event_year := year(peak_date)]
  
  if (save_event_panel) {
    out_path_event <- file.path(out_dir, paste0("panel_", region_tag, ".rds"))
    saveRDS(panel2, out_path_event)
    message("  Saved event panel → ", out_path_event)
  }
  
  invisible(list(panel = panel, panel2 = panel2))
}

# =============================================================
# Call
# =============================================================
build_panel(
  target_huc4s     = "1802",
  region_tag       = "sac",
  upstream_file    = "weighted_lc_runoff_1802.rds",
  T_pre            = 7L,
  T_post           = 7L,
  save_main_panel  = FALSE,
  save_event_panel = TRUE
)
