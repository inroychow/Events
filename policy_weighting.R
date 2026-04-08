library(data.table)
library(stringr)
library(tidyverse)

setDT(panel2)

# Load needed data in here ----------
huc12_landcover = readRDS("sac/data/lc_huc12_sac.rds") #load in landcover sac, just share developed and share undeveloped (not weighted by area rn)
upstream_map= readRDS("sac/data/htdt_sac.rds") #huc12 to downstream tract mapping
panel2 = readRDS("sac/data/panel2.rds")
setDT(panel2)
setDT(upstream_map)
setDT(huc12_landcover)
# ------------------------------------

# one row per tract-year for policy baseline
tract_2009_policies <- panel2[year == 2009,
                              .(active_policies_2009 = sum(active_policies, na.rm = TRUE)),
                              by = .(tract)]
tract_2009_policies[, tract := str_pad(as.character(tract), width = 11, side = "left", pad = "0")] #pad 


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

huc12_landcover_w <- merge(
  huc12_landcover,
  huc12_weights,
  by = "huc12",
  all.x = TRUE
)

huc12_landcover_w[, weight_2009 := fifelse(is.na(weight_2009), 0, weight_2009)]

huc12_landcover_w[, huc4:="1802"] #make row to indicate huc4

#Weight land cover by downstream policies
huc4_landcover <- huc12_landcover_w[
  ,
  .(
    weighted_developed_share =
      if (sum(weight_2009, na.rm = TRUE) > 0)
        sum(developed_share * weight_2009, na.rm = TRUE) / sum(weight_2009, na.rm = TRUE)
    else NA_real_,
    
    weighted_undeveloped_share =
      if (sum(weight_2009, na.rm = TRUE) > 0)
        sum(undeveloped_share * weight_2009, na.rm = TRUE) / sum(weight_2009, na.rm = TRUE)
    else NA_real_,
    
    weighted_developed_sqkm =
      if (sum(weight_2009, na.rm = TRUE) > 0)
        sum(developed_sqkm * weight_2009, na.rm = TRUE) / sum(weight_2009, na.rm = TRUE)
    else NA_real_,
    
    weighted_undeveloped_sqkm =
      if (sum(weight_2009, na.rm = TRUE) > 0)
        sum(undeveloped_sqkm * weight_2009, na.rm = TRUE) / sum(weight_2009, na.rm = TRUE)
    else NA_real_,
    
    total_weight_2009 = sum(weight_2009, na.rm = TRUE)
  ),
  by = .(huc4, year)
]

panel2[,huc4:="1802"]

# Aggregate panel2 to HUC4-date level 
# Sum paid and coverage across all tracts in the basin, per day
huc4_panel <- panel2[,
                     .(
                       total_paid       = sum(tract_total_paid,  na.rm = TRUE),
                       total_coverage   = sum(active_coverage,   na.rm = TRUE),
                       total_claims     = sum(tract_claims,      na.rm = TRUE),
                       total_policies   = sum(active_policies,   na.rm = TRUE),
                       
                       # event structure: all tracts share the same event (tau) within a day,
                       # so just take the first non-NA value
                       tau              = first(na.omit(tau)),
                       in_event_window  = first(na.omit(in_event_window)),
                       event_id         = first(na.omit(event_id)),
                       event_start      = first(na.omit(event_start)),
                       event_end        = first(na.omit(event_end)),
                       extreme_event    = first(na.omit(extreme_event)),
                       peak_date        = first(na.omit(peak_date)),
                       peak_precip      = first(na.omit(peak_precip)),
                       precip_mm        = mean(precip_mm,         na.rm = TRUE),  # or sum — your call
                       roll3            = mean(roll3,             na.rm = TRUE),
                       is_peak_day      = max(is_peak_day,        na.rm = TRUE)
                     ),
                     by = .(huc4, date, year, month, yday, dow, event_year)
]

# Outcome, Y_dy: basin-wide damages as fraction of total coverage ─────────────
huc4_panel[, Y_dy := fifelse(total_coverage > 0,
                             total_paid / total_coverage,
                             NA_real_)]

huc4_panel[, log_Y_dy := fifelse(Y_dy > 0, log(Y_dy), NA_real_)]

# Tau factor (for dummies W_tau) 
# Re-create tau_factor at huc4 level 
tau_levels <- as.character(-7:7)  
huc4_panel[, tau_factor := factor(tau, levels = tau_levels)]

# Merge in basin-wide weighted land cover (huc4 × year) 
# huc4_landcover already has: huc4, year, weighted_developed_share,
#   weighted_undeveloped_share, weighted_developed_sqkm, weighted_undeveloped_sqkm

huc4_panel <- merge(
  huc4_panel,
  huc4_landcover[, .(huc4, year,
                     weighted_developed_share,
                     weighted_undeveloped_share,
                     weighted_developed_sqkm,
                     weighted_undeveloped_sqkm)],
  by  = c("huc4", "year"),
  all.x = TRUE
)


head(huc4_panel)
