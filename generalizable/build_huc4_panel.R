library(data.table)
library(stringr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(fixest)
setwd("L:/Wetland Flood Mitigation/Paper_NFIP")

# -----------------------------
# Load needed data
# -----------------------------
huc12_landcover <- readRDS("sac/data/lc_huc12_sac.rds")
upstream_map    <- readRDS("sac/data/htdt_sac.rds")
panel2          <- readRDS("sac/data/panel2_sac.rds")
precip_huc12    <- readRDS("data/precip_huc12_tiles/precip_huc12s_all.rds")

precip_huc12 <- precip_huc12 %>%
  rename(huc12 = upstream_huc12)

setDT(precip_huc12)
setDT(panel2)
setDT(upstream_map)
setDT(huc12_landcover)

# -----------------------------
# Standardize IDs
# -----------------------------
if ("tract" %in% names(panel2)) {
  panel2[, tract := str_pad(as.character(tract), width = 11, side = "left", pad = "0")]
}
if ("tract" %in% names(upstream_map)) {
  upstream_map[, tract := str_pad(as.character(tract), width = 11, side = "left", pad = "0")]
}
if ("huc12" %in% names(huc12_landcover)) {
  huc12_landcover[, huc12 := str_pad(as.character(huc12), width = 12, side = "left", pad = "0")]
}
if ("upstream_huc12" %in% names(upstream_map)) {
  upstream_map[, upstream_huc12 := str_pad(as.character(upstream_huc12), width = 12, side = "left", pad = "0")]
}



# 2009 policy weights per tract
tract_2009_policies <- panel2[
  year == 2009,
  .(active_policies_2009 = sum(active_policies, na.rm = TRUE)),
  by = tract
]

# sum policies across all downstream tracts for each HUC12
huc12_weights <- merge(upstream_map, tract_2009_policies, by = "tract", all.x = TRUE)
huc12_weights[, active_policies_2009 := fifelse(is.na(active_policies_2009), 0, active_policies_2009)]
huc12_weights <- huc12_weights[
  , .(weight_2009 = sum(active_policies_2009, na.rm = TRUE)),
  by = upstream_huc12
]
setnames(huc12_weights, "upstream_huc12", "huc12")

# weighted precip at HUC4-day
precip_w <- merge(precip_huc12, huc12_weights, by = "huc12", all.x = TRUE)
precip_w[, weight_2009 := fifelse(is.na(weight_2009), 0, weight_2009)]
precip_w[, huc4 := substr(huc12, 1, 4)]
precip_huc4 <- precip_w[
  ,
  .(
    precip_mm_wtd = fifelse(
      sum(weight_2009) > 0,
      weighted.mean(precip_mm, w = weight_2009, na.rm = TRUE),
      mean(precip_mm, na.rm = TRUE)
    )
  ),
  by = .(huc4, date)
]

panel2[, huc4 := "1802"]

huc4_panel <- panel2[
  ,
  .(
    total_paid               = sum(tract_total_paid,         na.rm = TRUE),
    total_coverage           = sum(active_coverage,          na.rm = TRUE),
    total_claims             = sum(tract_claims,             na.rm = TRUE),
    total_policies           = sum(active_policies,          na.rm = TRUE),
    upstream_runoff_in       = mean(upstream_runoff_in,      na.rm = TRUE),
    upstream_precip_3day_in  = mean(upstream_precip_3day_in, na.rm = TRUE),
    developed_share          = mean(developed_share,         na.rm = TRUE),
    wetland_share            = mean(wetland_share,           na.rm = TRUE),
    forest_share             = mean(forest_share,            na.rm = TRUE),
    undeveloped_share        = mean(undeveloped_share,       na.rm = TRUE),
    tau             = first(na.omit(tau)),
    in_event_window = first(na.omit(in_event_window)),
    event_id        = first(na.omit(event_id)),
    event_start     = first(na.omit(event_start)),
    event_end       = first(na.omit(event_end)),
    extreme_event   = first(na.omit(extreme_event)),
    peak_date       = first(na.omit(peak_date)),
    peak_precip     = first(na.omit(peak_precip)),
    is_peak_day     = max(is_peak_day, na.rm = TRUE)
  ),
  by = .(huc4, date, year, month, yday, dow, event_year)
]

# merge weighted precip
huc4_panel <- merge(huc4_panel, precip_huc4, by = c("huc4", "date"), all.x = TRUE)

# outcomes
huc4_panel[, Y_dy     := fifelse(total_coverage > 0, total_paid / total_coverage, NA_real_)]
huc4_panel[, log_Y_dy := fifelse(Y_dy > 0, log(Y_dy), NA_real_)]

# tau factor
tau_levels <- as.character(-7:7)
huc4_panel[, tau_factor := relevel(factor(as.character(tau), levels = tau_levels), ref = "-1")]

