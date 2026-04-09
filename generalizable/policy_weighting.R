library(data.table)
library(stringr)
library(tidyverse)
library(data.table)
library(tidyr)
library(dplyr)
# Load needed data in here ----------
huc12_landcover = readRDS("sac/data/lc_huc12_sac.rds") #load in landcover sac, just share developed and share undeveloped (not yet area weighted)
upstream_map= readRDS("sac/data/htdt_sac.rds") #huc12 to downstream tract mapping
panel2 = readRDS("sac/data/panel2_sac.rds")
setDT(panel2)
setDT(upstream_map)
setDT(huc12_landcover)
# ------------------------------------

# start from lc_shares
huc12_landcover <- as.data.table(huc12_landcover)

# make long
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
  mutate(year = as.integer(year)) %>%
  as.data.table()

huc12_long_w <- merge(
  huc12_long,
  huc12_weights,   # columns: huc12, weight_2009
  by = "huc12",
  all.x = TRUE
)

# Merge baseline policies onto downstream tract map
huc12_weights <- merge(
  upstream_map,
  tract_2009_policies,
  by = "tract",
  all.x = TRUE
)

# Replace missing downstream policy counts with 0
huc12_weights[, active_policies_2009 := fifelse(is.na(active_policies_2009), 0, active_policies_2009)]

# Fixed HUC12 weight = total downstream baseline policies
huc12_weights <- huc12_weights[
  , .(weight_2009 = sum(active_policies_2009, na.rm = TRUE)),
  by = .(upstream_huc12)
]

huc12_weights = huc12_weights %>% 
  rename(huc12 = upstream_huc12)

huc12_long_w <- merge(
  huc12_long,
  huc12_weights,   # columns: huc12, weight_2009
  by = "huc12",
  all.x = TRUE
)

huc12_long_w[, weight_2009 := fifelse(is.na(weight_2009), 0, weight_2009)]


huc4_landcover <- huc12_long_w[
  ,
  .(
    weighted_developed_sqkm =
      sum(developed_sqkm * weight_2009, na.rm = TRUE),
    
    weighted_undeveloped_sqkm =
      sum(undeveloped_sqkm * weight_2009, na.rm = TRUE),
    
    weighted_total_sqkm =
      sum(area_sqkm_huc12 * weight_2009, na.rm = TRUE)
  ),
  by = .(huc4, year)
]

huc4_landcover[, weighted_developed_share :=
                 fifelse(weighted_total_sqkm > 0,
                         weighted_developed_sqkm / weighted_total_sqkm,
                         NA_real_)]

huc4_landcover[, weighted_undeveloped_share :=
                 fifelse(weighted_total_sqkm > 0,
                         weighted_undeveloped_sqkm / weighted_total_sqkm,
                         NA_real_)]

###

panel2[,huc4:="1802"]


# Aggregate panel2 to HUC4-date level
huc4_panel <- panel2[
  ,
  .(
    total_paid      = sum(tract_total_paid, na.rm = TRUE),
    total_coverage  = sum(active_coverage, na.rm = TRUE),
    total_claims    = sum(tract_claims, na.rm = TRUE),
    total_policies  = sum(active_policies, na.rm = TRUE),
    
    tau             = first(na.omit(tau)),
    in_event_window = first(na.omit(in_event_window)),
    event_id        = first(na.omit(event_id)),
    event_start     = first(na.omit(event_start)),
    event_end       = first(na.omit(event_end)),
    extreme_event   = first(na.omit(extreme_event)),
    peak_date       = first(na.omit(peak_date)),
    peak_precip     = first(na.omit(peak_precip)),
    precip_mm       = mean(precip_mm, na.rm = TRUE),
    roll3           = mean(roll3, na.rm = TRUE),
    is_peak_day     = max(is_peak_day, na.rm = TRUE)
  ),
  by = .(huc4, date, year, month, yday, dow, event_year)
]

# Outcome
huc4_panel[, Y_dy := fifelse(total_coverage > 0,
                             total_paid / total_coverage,
                             NA_real_)]

huc4_panel[, log_Y_dy := fifelse(Y_dy > 0, log(Y_dy), NA_real_)]

# Tau factor
tau_levels <- as.character(-7:7)
huc4_panel[, tau_factor := factor(as.character(tau), levels = tau_levels)]
huc4_panel[, tau_factor := relevel(tau_factor, ref = "-1")]

# Make sure year is integer
huc4_landcover[, year := as.integer(year)]

# Get all weighted land-cover columns automatically
weighted_cols <- grep("^weighted_", names(huc4_landcover), value = TRUE)

# Keep only huc4, year, and all weighted variables
nlcd_years <- unique(
  huc4_landcover[, c("huc4", "year", weighted_cols), with = FALSE]
)

# Set keys for rolling join
setkey(nlcd_years, huc4, year)
setkey(huc4_panel, huc4, year)

# Rolling join: assign each panel year the most recent NLCD year <= it
huc4_panel <- nlcd_years[huc4_panel, roll = TRUE]


saveRDS(huc4_panel, "sac/data/huc4_panel_sac.rds")

m1 <- feols(
  Y_dy ~ tau_factor + tau_factor : weighted_developed_share
  | event_id,
  data = huc4_panel[in_event_window == 1],
  vcov = "hetero"
)

etable(m1)
