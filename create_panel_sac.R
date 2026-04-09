library(data.table)
library(dplyr)

# Run build_lc_sac.R build_nfip_sac.R and build_events_sac_2009.R first.

claims_sac = readRDS("sac/data/claims_sac.rds")
policies_sac = readRDS("sac/data/policies_sac.rds")
event_days_post2009 = readRDS("sac/data/eventmap_post2009_sac.rds")
tract_year_lc = readRDS("data/tract_year_lc2.rds")
#-----------------------------#
# 1. Standardize inputs
#-----------------------------#

setDT(claims_sac)
setDT(policies_sac)
setDT(event_days_post2009)

tract_year_lc_dt <- as.data.table(tract_year_lc)

# Make sure IDs are character
claims_sac[, tract := as.character(tract)]
policies_sac[, tract := as.character(tract)]
tract_year_lc_dt[, tract := as.character(tract)]

# Convert dates
claims_sac[, dateOfLoss := as.Date(dateOfLoss)]
policies_sac[, policyEffectiveDate := as.Date(policyEffectiveDate)]
policies_sac[, policyTerminationDate := as.Date(policyTerminationDate)]

event_days_post2009[, event_start := as.Date(event_start)]
event_days_post2009[, event_end   := as.Date(event_end)]
event_days_post2009[, peak_date   := as.Date(peak_date)]

# land cover year as integer
tract_year_lc_dt[, year := as.integer(year)]

#-----------------------------#
# 2. Aggregate claims to tract-day
#-----------------------------#

claims_day <- claims_sac[
  , .(
    tract_claims = sum(tract_claims, na.rm = TRUE),
    tract_total_paid = sum(tract_total_paid, na.rm = TRUE)
  ),
  by = .(tract, date = dateOfLoss)
]

#-----------------------------#
# 3. Build full tract-day panel
#-----------------------------#

# tract universe: any tract appearing in claims, policies, or land cover
all_tracts <- unique(c(
  claims_sac$tract,
  policies_sac$tract,
  tract_year_lc_dt$tract
))

panel_start <- as.Date("2009-01-01")
# study dates:
# start at earliest relevant post-2009 date shared by your data
panel_end <- max(
  max(claims_sac$dateOfLoss, na.rm = TRUE),
  max(policies_sac$policyTerminationDate, na.rm = TRUE),
  max(event_days_post2009$event_end, na.rm = TRUE),
  na.rm = TRUE
)

panel_dates <- seq(panel_start, panel_end, by = "day")

panel <- CJ(
  tract = all_tracts,
  date = panel_dates,
  unique = TRUE
)

#-----------------------------#
# 4. Join claims onto panel
#-----------------------------#

panel <- claims_day[panel, on = .(tract, date)]

# Fill missing claims with zero
panel[is.na(tract_claims), tract_claims := 0L]
panel[is.na(tract_total_paid), tract_total_paid := 0]


#-----------------------------#
# 5. Add active policies by tract-day
#-----------------------------#
#-----------------------------#
# 5. Add active policies by tract-day
#-----------------------------#
policy_intervals <- policies_sac[
  , .(
    tract,
    start = policyEffectiveDate,
    end   = policyTerminationDate,
    total_coverage
    # drop n_policies_tract — it's a static tract total, not useful here
  )
]

policy_daily <- policy_intervals[
  panel,
  .(tract, date = i.date, total_coverage),
  on  = .(tract, start <= date, end >= date),
  allow.cartesian = TRUE,
  nomatch = 0L
][
  , .(
    active_policies = .N,                                 # policies active on this date
    active_coverage = sum(total_coverage, na.rm = TRUE)   # total insured value active
  ),
  by = .(tract, date)
]

panel <- policy_daily[panel, on = .(tract, date)]
panel[is.na(active_policies), active_policies := 0L]
panel[is.na(active_coverage), active_coverage := 0]

#-----------------------------#
# 6. Add extreme-event indicator
#-----------------------------#

event_daily <- unique(
  event_days_post2009[, .(
    date,
    extreme_event,
    event_id,
    event_start,
    event_end,
    peak_date,
    peak_precip,
    precip_mm,
    roll3,
    is_peak_day
  )]
)

panel <- event_daily[panel, on = "date"]

panel[is.na(extreme_event), extreme_event := 0L]
panel[is.na(is_peak_day), is_peak_day := 0L]

#-----------------------------#
# 7. Assign each tract-day a land cover year
#-----------------------------#

# Most recent available NLCD snapshot
panel[, lc_year :=
        fifelse(date <  as.Date("2012-01-01"), 2009L,
                fifelse(date <  as.Date("2015-01-01"), 2012L,
                        fifelse(date <  as.Date("2018-01-01"), 2015L,
                                fifelse(date <  as.Date("2021-01-01"), 2018L,
                                        fifelse(date <  as.Date("2024-01-01"), 2021L, 2024L)))))]


#-----------------------------#
# 8. Join land cover to tract-day panel
#-----------------------------#

# rename year to lc_year for join
setnames(tract_year_lc_dt, "year", "lc_year")

panel <- tract_year_lc_dt[
  panel,
  on = .(tract, lc_year)
]

# panel is now enriched with:
# - claims
# - active policies
# - event indicator
# - land cover measures
#-----------------------------#
# 9. Add useful time variables and outcomes
#-----------------------------#

panel[
  , `:=`(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    yday = as.integer(format(date, "%j")),
    dow = weekdays(date)
  )
]

# Common transformed outcomes
panel[
  , `:=`(
    log_total_paid = log1p(tract_total_paid),
    log_active_average = log1p(active_coverage)
  )
]

#-----------------------------#
# 10. Optional: create interaction terms now
#-----------------------------#

landcover_vars <- c(
  "agriculture_share", "barren_share", "developed_share",
  "forest_share", "shrub_grass_share", "water_share", "wetland_share"
)

for (v in landcover_vars) {
  new_name <- paste0("event_x_", v)
  panel[, (new_name) := get(v) * extreme_event]
}

panel[, claims_per_coverage := fifelse(active_coverage > 0, tract_total_paid / active_coverage, NA_real_)]
saveRDS(panel, "sac/data/panel_sac2.rds")
saveRDS(panel, "data/panel_sac_all.rds")
