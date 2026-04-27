setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
# Load packages
pacman::p_load(sf, terra, exactextractr, data.table, dplyr, lubridate, stringr, prism, fixest, future, future.apply, progressr, tidyverse, fst, patchwork, scales)

# ----------------- THREE DAY THRESHOLD

# First, run create_nfip_lc_sac.R and extract_precip_sac.R
nfip_htdt_sac_lc = readRDS("data/lcchange_nfip_sac.rds")

precip_all = readRDS("data/precip_huc12_tiles/precip_huc12s_all.rds")

precip_daily <- precip_all[, .(precip_mm = mean(precip_mm, na.rm = TRUE)), by = date] # use mean or max?

# =============================================================================
# PRECIP → NFIP EVENT LINKAGE — SACRAMENTO BASIN
# =============================================================================


# ── 1. Rolling 3-day total ─────────────────────────────────────────────────
precip_daily[order(date),
             roll3 := rollsum(precip_mm, k = 3, fill = NA, align = "right")]

# ── 2. 99th-pct threshold on RAINY days only ──────────────────────────────
threshold <- quantile(
  precip_daily[precip_mm > 1, roll3],
  0.99, na.rm = TRUE
)
cat(sprintf("Event threshold (99th pct, rainy days): %.1f mm over 3 days\n", threshold))

# ── 3. Flag raw event dates ────────────────────────────────────────────────
precip_daily[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]

# ── 4. Expand ±3 day window around each event date ────────────────────────
#   For each flagged date, mark date-3 through date+3 as "in window"
event_dates <- precip_daily[precip_event_raw == TRUE, date]

window_dates <- unique(as.Date(unlist(lapply(event_dates, function(d) {
  seq(d - 3, d + 3, by = "day")
}))))

precip_daily[, in_event_window := date %in% window_dates]

# ── 5. Collapse contiguous windows into discrete EVENT IDs ────────────────
#   Avoids double-counting when consecutive days exceed threshold
precip_daily[order(date),
             event_block := cumsum(c(FALSE, diff(in_event_window * 1) == 1))]

# Only keep event blocks where in_event_window is TRUE
precip_daily[in_event_window == FALSE, event_block := NA]

# Assign clean event IDs (1, 2, 3...)
event_map <- precip_daily[!is.na(event_block),
                          .(event_start = min(date),
                            event_end   = max(date),
                            peak_precip = max(roll3, na.rm = TRUE),
                            peak_date   = date[which.max(roll3)]),
                          by = event_block][
                            , event_id := .I][]

# IDENTIFY PRECIP EVENTS THAT EXCEED THE 99th percentile THRESHOLD---------

cat(sprintf("Identified %d discrete precip events\n", nrow(event_map)))

saveRDS(event_map, "data/eventmap.rds")

# ---------------- past 2000 only 

precip_post2000 <- precip_daily[date >= as.Date("2000-01-01")]

threshold <- quantile(
  precip_post2000[precip_mm > 1, roll3],
  0.99,
  na.rm = TRUE
)

cat(sprintf("Event threshold (99th pct, rainy days, 2000+): %.1f mm over 3 days\n", threshold))

precip_daily[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]

# subset to post-2000 only
precip_post2000 <- copy(precip_daily[date >= as.Date("2000-01-01")])

# flag events using threshold defined from post-2000 data
precip_post2000[, precip_event_raw := (!is.na(roll3)) & (roll3 >= threshold)]

# dates that exceed threshold
event_dates <- precip_post2000[precip_event_raw == TRUE, date]

# expand each event by +/- 3 days
window_dates <- unique(as.Date(unlist(lapply(event_dates, function(d) {
  seq(d - 3, d + 3, by = "day")
}))))

# mark whether each day is inside an event window
precip_post2000[, in_event_window := date %in% window_dates]

# assign blocks of contiguous event-window days
precip_post2000[order(date),
                event_block := cumsum(c(FALSE, diff(in_event_window * 1) == 1))]

# set non-event days to NA
precip_post2000[in_event_window == FALSE, event_block := NA]

# create event map
event_map_post2000 <- precip_post2000[!is.na(event_block),
                                      .(event_start = min(date),
                                        event_end   = max(date),
                                        peak_precip = max(roll3, na.rm = TRUE),
                                        peak_date   = date[which.max(roll3)]),
                                      by = event_block][, event_id := .I][]

event_map_post2000

saveRDS(event_map_post2000, "data/eventmap_post2000.rds")

####################################################################


# ── 6. Assign each calendar date to an event_id (NA = non-event) ──────────
#   Using a non-equi join: date falls within [event_start, event_end]
precip_daily[, event_id := NA_integer_]
for (e in seq_len(nrow(event_map))) {
  precip_daily[date %between% c(event_map$event_start[e], event_map$event_end[e]),
               event_id := event_map$event_id[e]]
}

# ── 7. Join event flag onto NFIP panel ────────────────────────────────────
#nfip_htdt_sac_lc is currently later in the code ; need to MOVE IT
nfip_htdt_sac_lc[, loss_date := as.Date(dateOfLoss)]
nfip_event <- precip_daily[, .(date, roll3, in_event_window, event_id)][
  nfip_htdt_sac_lc[!is.na(loss_date)],   # drop the NA loss_date rows
  on = c("date" = "loss_date")
]

nfip_event[, in_event := !is.na(event_id)]


# ── 8. Summarise: event vs non-event damage ────────────────────────────────
nfip_event[, .(
  n_claims     = sum(tract_claims,      na.rm = TRUE),
  total_damage = sum(tract_total_paid,  na.rm = TRUE),
  n_days       = uniqueN(dateOfLoss)
), by = in_event][
  , avg_daily_damage := total_damage / n_days][]

# ── 9. Event-level summary: damage by event ───────────────────────────────
event_damage <- nfip_event[in_event == TRUE, .(
  n_claims     = sum(tract_claims,     na.rm = TRUE),
  total_damage = sum(tract_total_paid, na.rm = TRUE)
), by = event_id][
  event_map, on = "event_id"
][order(event_start)]

print(event_damage[, .(event_id, event_start, event_end,
                       peak_precip, n_claims, total_damage)])

# ------------------------------------

# ── Prep ─────────────────────────────────────────────────────────────────────
ed <- copy(event_damage)
ed[, year := year(event_start)]
ed[, log_damage := log1p(total_damage)]
ed[is.na(total_damage), total_damage := 0]
ed[is.na(n_claims),     n_claims     := 0]

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_flood <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(color = "grey50", size = 9),
    axis.title         = element_text(size = 9)
  )

# ── 1. Timeline: damage per event, sized by claims ───────────────────────────
p_precipevents <- ggplot(ed[total_damage > 0],
                         aes(x = event_start, y = total_damage,
                             size = n_claims, color = peak_precip)) +
  geom_segment(aes(xend = event_start, yend = 0),
               linewidth = 0.3, color = "grey70") +
  geom_point(alpha = 0.75) +
  scale_y_continuous(labels = dollar_format(scale = 1e-6, suffix = "M"),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_color_distiller(palette = "Blues", direction = 1,
                        name = "Peak 3-day\nprecip (mm)") +
  scale_size_continuous(range = c(1, 8), name = "# claims",
                        labels = comma) +
  labs(title    = "NFIP Damage per Precip Event — Sacramento Basin",
       x = NULL, y = "Total damage paid") +
  theme_flood
p_precipevents

# saveRDS(p1, "Figures/precipevents_sac.png")


# TOP 5 EVENTS ----------------------------------------------------------------------

# ── Pick top 5 events by total damage ─────────────────────────────────────
top5 <- event_damage[order(-total_damage)][1:5] #all data
top5 <- event_damage[year(event_start) >= 2000][order(-total_damage)][1:5] # just past 2000

# ── Build long data for all 5 events ──────────────────────────────────────
precip_long <- rbindlist(lapply(seq_len(nrow(top5)), function(i) {
  w_start <- top5$event_start[i] - 7
  w_end   <- top5$event_end[i]   + 7
  dt <- precip_daily[date %between% c(w_start, w_end)]
  dt[, `:=`(
    event_label = sprintf("%s – %s",
                          format(top5$event_start[i], "%b %d %Y"),
                          format(top5$event_end[i],   "%b %d %Y")),
    event_start = top5$event_start[i],
    event_end   = top5$event_end[i],
    peak_date   = top5$peak_date[i],
    peak_precip = top5$peak_precip[i],
    day         = as.integer(date - top5$event_start[i])
  )]
}))

nfip_long <- rbindlist(lapply(seq_len(nrow(top5)), function(i) {
  eid     <- top5$event_id[i]
  w_start <- top5$event_start[i] - 7
  w_end   <- top5$event_end[i]   + 7
  dt <- nfip_event[!is.na(date) &
                     date %between% c(w_start, w_end) &
                     event_id == eid,
                   .(daily_damage = sum(tract_total_paid, na.rm = TRUE),
                     daily_claims = sum(tract_claims,     na.rm = TRUE)),
                   by = date]
  dt[, `:=`(
    event_label = sprintf("%s – %s",
                          format(top5$event_start[i], "%b %d %Y"),
                          format(top5$event_end[i],   "%b %d %Y")),
    event_start = top5$event_start[i],
    event_end   = top5$event_end[i],
    day         = as.integer(date - top5$event_start[i])
  )]
}))

# ── Rect data in relative days ─────────────────────────────────────────────
rect_long <- unique(precip_long[, .(
  event_label,
  event_start_day = 0L,
  event_end_day   = as.integer(event_end - event_start)
)])

# ── Peak precip day (relative) ─────────────────────────────────────────────
peak_long <- unique(precip_long[, .(
  event_label,
  peak_day = as.integer(peak_date - event_start)
)])

precip_long <- peak_long[precip_long, on = "event_label"]

common_x <- scale_x_continuous(
  limits = c(-7, 30),
  breaks = seq(-7, 30, by = 7),
  labels = function(x) paste0("d", x),
  expand = c(0, 0)
)

p_top5_precip <- ggplot(precip_long, aes(x = day, y = precip_mm)) +
  geom_rect(data = rect_long,
            aes(xmin = event_start_day, xmax = event_end_day,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "steelblue", alpha = 0.12) +
  geom_area(fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue") +
  geom_vline(aes(xintercept = peak_day),
             linetype = "dashed", color = "navy", linewidth = 0.4) +
  facet_wrap(~ event_label, nrow = 1) +
  common_x +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Daily precip (mm)",
       title = "Top 5 Precip Events — Sacramento Basin (2000 onwards)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11))

p_top5_claims <- ggplot(nfip_long, aes(x = day, y = daily_damage)) +
  geom_rect(data = rect_long,
            aes(xmin = event_start_day, xmax = event_end_day,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "coral", alpha = 0.12) +
  geom_col(fill = "coral", alpha = 0.8, width = 1) +
  facet_wrap(~ event_label, nrow = 1) +
  common_x +
  scale_y_continuous(labels = dollar_format(scale = 1e-3, suffix = "K")) +
  labs(x = "Days since event start", y = "Damage paid",
       title = "NFIP Claims During Events") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11))

p_top5_precip / p_top5_claims

# ── Save ───────────────────────────────────────────────────────────────────
# ggsave("Figures/sactop5.png", width = 16, height = 7)
# ggsave("Figures/sactop5_past2000.png", width = 16, height = 7)