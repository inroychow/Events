setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
# Load packages
pacman::p_load(sf, terra, exactextractr, data.table, dplyr, lubridate, stringr, prism, fixest, future, future.apply, progressr, tidyverse, fst, patchwork, scales)

# ----------------- THREE DAY THRESHOLD

# First, run create_nfip_lc_sac.R and extract_precip_sac.R

precip_all = readRDS("data/precip_huc12_tiles/precip_huc12s_all.rds")

precip_daily <- precip_all[, .(precip_mm = mean(precip_mm, na.rm = TRUE)), by = date] # use mean or max?

# =============================================================================
# PRECIP → NFIP EVENT LINKAGE — SACRAMENTO BASIN
# =============================================================================


# ── 1. Rolling 3-day total ─────────────────────────────────────────────────
precip_daily[order(date),
             roll3 := rollsum(precip_mm, k = 3, fill = NA, align = "right")]

# ---------------- past 2009 only 

precip_post2009 <- precip_daily[date >= as.Date("2009-01-01")]

threshold <- quantile(
  precip_post2009[precip_mm > 1, roll3],
  0.99,
  na.rm = TRUE
) # 99th percentile events

cat(sprintf("Event threshold (99th pct, rainy days, 2009+): %.1f mm over 3 days\n", threshold))

precip_daily[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]

# subset to post-2009 only
precip_post2009 <- copy(precip_daily[date >= as.Date("2009-01-01")])

# flag events using threshold defined from post-2009 data
precip_post2009[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]

# dates that exceed threshold
event_dates <- precip_post2009[precip_event_raw == TRUE, date]

# expand each event by +/- 3 days
window_dates <- unique(as.Date(unlist(lapply(event_dates, function(d) {
  seq(d - 3, d + 3, by = "day")
}))))

# mark whether each day is inside an event window
precip_post2009[, in_event_window := date %in% window_dates]

# assign blocks of contiguous event-window days
precip_post2009[order(date),
                event_block := cumsum(c(FALSE, diff(in_event_window * 1) == 1))]

# set non-event days to NA
precip_post2009[in_event_window == FALSE, event_block := NA]

# create event map
event_map_post2009 <- precip_post2009[!is.na(event_block),
                                      .(event_start = min(date),
                                        event_end   = max(date),
                                        peak_precip = max(roll3, na.rm = TRUE),
                                        peak_date   = date[which.max(roll3)]),
                                      by = event_block][, event_id := .I][]

event_map_post2009

event_map_post2009 <- precip_post2009[!is.na(event_block),
                                      .(event_start = min(date),
                                        event_end   = max(date),
                                        peak_precip = max(roll3, na.rm = TRUE),
                                        peak_date   = date[which.max(roll3)]),
                                      by = event_block][, event_id := .I][]

saveRDS(event_map_post2009, "data/eventmap_post2009.rds")

###

# Make daily event table: one row per day in each event window
# Requires: precip_post2009 and event_map_post2009 already exist

library(data.table)

setDT(precip_post2009)
setDT(event_map_post2009)

event_days_post2009 <- precip_post2009[!is.na(event_block),
                                       .(date, event_block, precip_mm, roll3)
][
  event_map_post2009[, .(
    event_block, event_id, event_start, event_end, peak_precip, peak_date
  )],
  on = "event_block"
][
  , `:=`(
    extreme_event = 1L,
    is_peak_day   = as.integer(date == peak_date)
  )
][
  order(event_id, date)
]

# optional save
saveRDS(event_days_post2009, "data/eventdays_post2009.rds")



