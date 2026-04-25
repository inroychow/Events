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
precip_huc12 = precip_huc12 %>% 
  rename(huc12=upstream_huc12)
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

# -----------------------------
# Build baseline tract weights from panel2
# using 2009 active policies
# -----------------------------
tract_2009_policies <- panel2[
  year == 2009,
  .(active_policies_2009 = sum(active_policies, na.rm = TRUE)),
  by = .(tract)
]

# -----------------------------
# Map tract baseline policies to upstream HUC12s
# -----------------------------
huc12_weights <- merge(
  upstream_map,
  tract_2009_policies,
  by = "tract",
  all.x = TRUE
)

huc12_weights[, active_policies_2009 := fifelse(is.na(active_policies_2009), 0, active_policies_2009)]

huc12_weights <- huc12_weights[
  ,
  .(weight_2009 = sum(active_policies_2009, na.rm = TRUE)),
  by = .(upstream_huc12)
]

setnames(huc12_weights, "upstream_huc12", "huc12")
huc12_weights[, huc12 := str_pad(as.character(huc12), width = 12, side = "left", pad = "0")]

# -----------------------------
# Build policy+area weighted precip at HUC4-day
# -----------------------------
precip_huc12[, huc12 := str_pad(as.character(huc12), width = 12, side = "left", pad = "0")]
precip_huc12[, huc4  := substr(huc12, 1, 4)]

huc12_area <- unique(huc12_landcover[, .(huc12, area_sqkm_huc12)])

precip_w <- merge(precip_huc12, huc12_weights, by = "huc12", all.x = TRUE)
precip_w <- merge(precip_w,     huc12_area,    by = "huc12", all.x = TRUE)

precip_w[, weight_2009     := fifelse(is.na(weight_2009),      0, weight_2009)]
precip_w[, area_sqkm_huc12 := fifelse(is.na(area_sqkm_huc12), 0, area_sqkm_huc12)]
precip_w[, combined_weight := weight_2009 * area_sqkm_huc12]

precip_huc4 <- precip_w[
  ,
  .(
    precip_mm_wtd = fifelse(
      sum(combined_weight, na.rm = TRUE) > 0,
      sum(precip_mm * combined_weight, na.rm = TRUE) / sum(combined_weight, na.rm = TRUE),
      mean(precip_mm, na.rm = TRUE)
    )
  ),
  by = .(huc4, date)
]

# -----------------------------
# Reshape HUC12 land cover to long by year
# -----------------------------
huc12_long <- huc12_landcover %>%
  pivot_longer(
    cols = matches("_(share|sqkm)_\\d{4}$"),
    names_to = c("landcover", "metric", "year"),
    names_pattern = "(.+)_(share|sqkm)_(\\d{4})",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols = c(huc12, huc4, area_sqkm_huc12, year),
    names_from = c(landcover, metric),
    values_from = value
  ) %>%
  mutate(
    year = as.integer(year),
    huc12 = str_pad(as.character(huc12), width = 12, side = "left", pad = "0"),
    huc4 = as.character(huc4)
  ) %>%
  as.data.table()

# -----------------------------
# Attach 2009 policy weights to HUC12 land cover
# -----------------------------
huc12_long_w <- merge(
  huc12_long,
  huc12_weights,
  by = "huc12",
  all.x = TRUE
)

huc12_long_w[, weight_2009 := fifelse(is.na(weight_2009), 0, weight_2009)]

# -----------------------------
# Aggregate to HUC4-year land cover
# policy-weighted and area-corrected
# -----------------------------
huc4_landcover <- huc12_long_w[
  ,
  .(
    weighted_developed_share =
      fifelse(
        sum(weight_2009 * area_sqkm_huc12, na.rm = TRUE) > 0,
        sum(developed_share * weight_2009 * area_sqkm_huc12, na.rm = TRUE) /
          sum(weight_2009 * area_sqkm_huc12, na.rm = TRUE),
        NA_real_
      ),
    weighted_undeveloped_share =
      fifelse(
        sum(weight_2009 * area_sqkm_huc12, na.rm = TRUE) > 0,
        sum(undeveloped_share * weight_2009 * area_sqkm_huc12, na.rm = TRUE) /
          sum(weight_2009 * area_sqkm_huc12, na.rm = TRUE),
        NA_real_
      )
  ),
  by = .(huc4, year)
]
# -----------------------------
# Build HUC4-day panel from panel2
# -----------------------------
panel2[, huc4 := "1802"]

huc4_panel <- panel2[
  ,
  .(
    total_paid      = sum(tract_total_paid,  na.rm = TRUE),
    total_coverage  = sum(active_coverage,   na.rm = TRUE),
    total_claims    = sum(tract_claims,      na.rm = TRUE),
    total_policies  = sum(active_policies,   na.rm = TRUE),
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

# -----------------------------
# Swap in policy+area weighted precip
# -----------------------------
huc4_panel <- merge(huc4_panel, precip_huc4, by = c("huc4", "date"), all.x = TRUE)

# -----------------------------
# Outcome variables
# -----------------------------
huc4_panel[, Y_dy := fifelse(total_coverage > 0,
                             total_paid / total_coverage,
                             NA_real_)]

huc4_panel[, log_Y_dy := fifelse(Y_dy > 0, log(Y_dy), NA_real_)]

# -----------------------------
# Tau factor
# -----------------------------
tau_levels <- as.character(-7:7)
huc4_panel[, tau_factor := factor(as.character(tau), levels = tau_levels)]
huc4_panel[, tau_factor := relevel(tau_factor, ref = "-1")]

# -----------------------------
# Join HUC4-year land cover onto daily panel
# using latest NLCD year <= panel year
# -----------------------------
huc4_landcover[, year := as.integer(year)]

weighted_cols <- grep("^weighted_", names(huc4_landcover), value = TRUE)

nlcd_years <- unique(
  huc4_landcover[, c("huc4", "year", weighted_cols), with = FALSE]
)

setkey(nlcd_years, huc4, year)
setkey(huc4_panel, huc4, year)

huc4_panel <- nlcd_years[huc4_panel, roll = TRUE]

#saveRDS(huc4_panel, "sac/data/huc4_panel_sac_wtd.rds")

# -----------------------------
# Model
# -----------------------------
m1 <- feols(
  Y_dy ~ tau_factor + tau_factor:weighted_developed_share
  | event_id,
  data = huc4_panel[in_event_window == 1],
  vcov = "hetero"
)

etable(m1)
