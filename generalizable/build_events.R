setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
# Load packages
pacman::p_load(sf, terra, exactextractr, data.table, dplyr, lubridate, stringr, prism, fixest, future, future.apply, progressr, tidyverse, fst, patchwork, scales)

# ----------------- THREE DAY THRESHOLD
# target_huc4s <- c("1806", "1801", "1805", "1804", "1803")
target_huc4s <- "1802"
# First, run extract_precip_sac.R

precip_all = readRDS("data/precip_huc12_tiles_ca/precip_huc12s_cahucs.rds")
setDT(precip_all)

# If huc4 is not already present, derive it from huc12
if (!"huc4" %in% names(precip_all)) {
  if ("huc12" %in% names(precip_all)) {
    precip_all[, huc4 := substr(as.character(huc12), 1, 4)]
  } else {
    stop("precip_all must contain either a 'huc4' column or a 'huc12' column.")
  }
}

# Keep only target HUC4s
precip_all <- precip_all[huc4 %in% target_huc4s]

# ------------------------------------------------------------
# FUNCTION TO BUILD EVENT TABLES FOR ONE HUC4
# ------------------------------------------------------------
build_event_tables_huc4 <- function(dt, huc4_id,
                                    start_date = as.Date("2009-01-01"),
                                    roll_k = 3,
                                    event_quantile = 0.97,
                                    rainy_day_cutoff = 1,
                                    pad_days = 3) {
  
  dt_huc4 <- copy(dt[huc4 == huc4_id])
  
  if (nrow(dt_huc4) == 0) {
    message("No data for HUC4 ", huc4_id)
    return(NULL)
  }
  
  # Collapse to daily precipitation for this HUC4
  precip_daily <- dt_huc4[, .(
    precip_mm = mean(precip_mm, na.rm = TRUE)
  ), by = date][order(date)]
  
  # Rolling 3-day precipitation
  precip_daily[, roll3 := zoo::rollsum(precip_mm, k = roll_k, fill = NA, align = "right")]
  
  # Post-2009 subset for threshold and event definition
  precip_post <- copy(precip_daily[date >= start_date])
  
  # Threshold based on rainy days only
  threshold <- quantile(
    precip_post[precip_mm > rainy_day_cutoff, roll3],
    probs = event_quantile,
    na.rm = TRUE
  )
  
  cat(sprintf(
    "HUC4 %s threshold (%.0fth pct, rainy days, %s+): %.1f mm over %d days\n",
    huc4_id, event_quantile * 100, format(start_date, "%Y"), threshold, roll_k
  ))
  
  # Flag threshold exceedance days
  precip_post[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]
  
  event_dates <- precip_post[precip_event_raw == TRUE, date]
  
  if (length(event_dates) == 0) {
    message("No events detected for HUC4 ", huc4_id)
    
    return(list(
      huc4 = huc4_id,
      threshold = threshold,
      precip_daily = precip_daily,
      precip_post = precip_post,
      event_map = data.table(),
      event_days = data.table()
    ))
  }
  
  # Expand each exceedance date by +/- pad_days
  window_dates <- unique(as.Date(unlist(lapply(event_dates, function(d) {
    seq(d - pad_days, d + pad_days, by = "day")
  }))))
  
  precip_post[, in_event_window := date %in% window_dates]
  
  # Assign contiguous event blocks
  precip_post[order(date),
              event_block := cumsum(c(FALSE, diff(in_event_window * 1) == 1))]
  
  precip_post[in_event_window == FALSE, event_block := NA_integer_]
  
  # Event summary table
  event_map <- precip_post[!is.na(event_block),
                           .(
                             event_start = min(date),
                             event_end   = max(date),
                             peak_precip = max(roll3, na.rm = TRUE),
                             peak_date   = date[which.max(roll3)]
                           ),
                           by = event_block][
                             , event_id := .I
                           ][
                             , huc4 := huc4_id
                           ][]
  
  # Daily event table
  event_days <- precip_post[!is.na(event_block),
                            .(date, event_block, precip_mm, roll3)][
                              event_map[, .(
                                event_block, event_id, event_start, event_end,
                                peak_precip, peak_date, huc4
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
  
  return(list(
    huc4 = huc4_id,
    threshold = threshold,
    precip_daily = precip_daily,
    precip_post = precip_post,
    event_map = event_map,
    event_days = event_days
  ))
}

# ------------------------------------------------------------
# RUN FOR ALL TARGET HUC4s
# ------------------------------------------------------------
results_by_huc4 <- lapply(target_huc4s, function(x) {
  build_event_tables_huc4(
    dt = precip_all,
    huc4_id = x,
    start_date = as.Date("2009-01-01"),
    roll_k = 3,
    event_quantile = 0.97,   # change to 0.95 or 0.99 if needed
    rainy_day_cutoff = 1,
    pad_days = 3
  )
})

names(results_by_huc4) <- target_huc4s

event_days_sac <- rbindlist(
  lapply(results_by_huc4, function(x) x$event_days),
  fill = TRUE
)

event_map_sac <- rbindlist(
  lapply(results_by_huc4, function(x) x$event_map),
  fill = TRUE
)

thresholds_by_huc4 <- rbindlist(lapply(results_by_huc4, function(x) {
  data.table(
    huc4 = x$huc4,
    threshold = x$threshold
  )
}), fill = TRUE)

# saveRDS(event_days_sac, "sac/data/eventdays_sac.rds")
