library(data.table)
library(fixest)
library(lubridate)
library(tidyverse)

panel = readRDS("sac/data/panel_sac.rds")

setDT(panel)
panel[, date := as.Date(date)]
panel[, event_start := as.Date(event_start)]
panel[, event_end := as.Date(event_end)]

# -----------------------------
# Keep unique events
# -----------------------------
events <- panel[!is.na(event_id) & !is.na(event_start),
             .(event_start = unique(event_start)[1],
               event_end   = unique(event_end)[1]),
             by = event_id]


# -----------------------------
# Choose event window... 7 days before/after event start date?
# -----------------------------
T_pre  <- 7
T_post <- 7

# expand each event into dates from event_start - T_pre to event_start + T_post
event_windows <- events[, .(
  date = seq(event_start - T_pre, event_start + T_post, by = "day")
), by = .(event_id, event_start, event_end)]

# relative day
event_windows[, tau := as.integer(date - event_start)]

# -----------------------------
# Attach event windows to all tracts
# -----------------------------
tracts <- unique(panel[, .(tract)])

event_windows[, join_id := 1L]
tracts[, join_id := 1L]

event_windows_tract <- merge(
  event_windows,
  tracts,
  by = "join_id",
  allow.cartesian = TRUE
)

event_windows_tract[, join_id := NULL]


# -----------------------------
# Merge event-window rows back to panel
# -----------------------------
panel2 <- merge(
  event_windows_tract,
  panel,
  by = c("tract", "date"),
  all.x = TRUE,
  allow.cartesian = TRUE
)

panel2 <- panel2[!is.na(date)]

# deal with name duplicate stuff
setnames(panel2,
         old = c("event_id.x", "event_start.x", "event_end.x"),
         new = c("event_id",   "event_start",   "event_end"))

# drop the original panel event-label variables
drop_cols <- intersect(c("event_id.y", "event_start.y", "event_end.y"), names(panel2))
panel2[, (drop_cols) := NULL]

# -----------------------------
# Create event-time indicators
# -----------------------------
panel2[, tau_factor := factor(tau)]

# omit tau = -1 as reference period (day before the event?)
panel2[, tau_factor := relevel(tau_factor, ref = "-1")]

# dummy for whether dates is inside actual event window
panel2[, in_event_window := fifelse(date >= event_start & date <= event_end, 1L, 0L)]
panel2[, event_year := year(event_start)]


# ----------------


library(fixest)

model <- feols(
  claims_per_coverage ~ i(tau, ref = -1) +                  # W_tau main effect
    i(tau, developed_share, ref = -1) +                     # W_tau × developed
    i(tau, undeveloped_share, ref = -1)                   # W_tau × undeveloped
   | tract,                                              # tract FE
  data   = panel2,
  cluster = ~tract
)

summary(model)
