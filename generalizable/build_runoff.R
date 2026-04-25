# =============================================================================
# Simplified Curve Number Construction for Sacramento HUC4 Basin
# =============================================================================
#
# Purpose: Build a simplified CN raster per NLCD panel year using only NLCD
#          land cover and one assumed hydrologic soil group, then area-average
#          to HUC12 subbasins.
#
# Main simplification:
#   - No gSSURGO / HSG raster
#   - CN is assigned directly from NLCD class using assumed HSG = "C"
#
# Outputs:
#   - data/cn/cn_raster_{year}.tif
#   - data/cn/cn_huc12.csv
#   - data/cn/runoff_huc12.csv
#   - data/upstream_runoff_tract.csv
#
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(fs)
library(slider)
library(dplyr)
library(tidyr)
library(readr)


setwd("L:/Wetland Flood Mitigation/Paper_NFIP")

# -----------------------------------------------------------------------------
# 1. FILE PATHS
# -----------------------------------------------------------------------------

data_dir <- "data"

nlcd_years <- c(2009, 2012, 2015, 2018, 2021, 2024)

nlcd_files <- setNames(
  file.path("L:/Wetland Flood Mitigation/Paper_NFIP/NLCD rasters", paste0("nlcd", nlcd_years, ".tif")),
  nlcd_years
)

cn_dir <- file.path(data_dir, "cn")
dir_create(cn_dir)

prism_dir <- file.path(data_dir, "prism", "daily")

# -----------------------------------------------------------------------------
# 2. LOAD AND PREPARE HUC12 BOUNDARIES
# -----------------------------------------------------------------------------

watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

target_huc4s <- "1802"

huc12_subset <- watersheds %>%
  mutate(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    huc4  = str_sub(huc12, 1, 4)
  ) %>%
  filter(huc4 %in% target_huc4s)

# Handle list columns if present
geom <- st_geometry(huc12_subset)

huc12_subset <- huc12_subset %>%
  st_drop_geometry() %>%
  mutate(across(
    where(is.list),
    ~ vapply(., function(x) paste(unlist(x), collapse = "|"), character(1))
  ))

st_geometry(huc12_subset) <- geom

# Use EPSG:5070 for raster work
huc12_sf <- huc12_subset %>%
  st_as_sf() %>%
  st_zm(drop = TRUE) %>%
  st_transform(5070) %>%
  mutate(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    huc4  = str_sub(huc12, 1, 4)
  )

# One basin boundary for cropping/masking
sac_huc4_buf <- huc12_sf %>%
  st_union() %>%
  st_as_sf()

cat(sprintf(
  "Selected HUC12s: %d polygons across %d HUC4s\n",
  nrow(huc12_sf),
  length(unique(huc12_sf$huc4))
))

# -----------------------------------------------------------------------------
# 3. CN LOOKUP TABLE
# -----------------------------------------------------------------------------

# cn_lookup <- tribble(
#   ~nlcd_code, ~nlcd_label,                      ~A,  ~B,  ~C,  ~D,
#   11L,        "Open water",                      98,  98,  98,  98,
#   21L,        "Developed, open space",            49,  69,  79,  84,
#   22L,        "Developed, low intensity",         57,  72,  81,  86,
#   23L,        "Developed, medium intensity",      61,  75,  83,  87,
#   24L,        "Developed, high intensity",        81,  88,  91,  93,
#   31L,        "Barren land",                      78,  86,  91,  93,
#   41L,        "Deciduous forest",                 45,  66,  77,  83,
#   42L,        "Evergreen forest",                 25,  55,  70,  77,
#   43L,        "Mixed forest",                     36,  60,  73,  79,
#   52L,        "Shrub/scrub",                      55,  72,  81,  86,
#   71L,        "Grassland/herbaceous",             50,  69,  79,  84,
#   81L,        "Pasture/hay",                      49,  69,  79,  84,
#   82L,        "Cultivated crops",                 67,  78,  85,  89,
#   90L,        "Woody wetlands",                   30,  58,  71,  78,
#   95L,        "Emergent herbaceous wetlands",     30,  58,  71,  78
# ) # https://www.hec.usace.army.mil/confluence/hmsdocs/hmsguides/gis-tools-and-terrain-data/gis-tutorials-and-guides/creating-a-curve-number-grid-and-computing-subbasin-average-curve-number-values


cn_lookup <- tribble(
  ~nlcd_code, ~nlcd_group,        ~A,  ~B,  ~C,  ~D,
  
  # Open water
  11L,        "Water",             98,  98,  98,  98,
  
  # Developed merged: average of 21, 22, 23, 24
  21L,        "Developed",         62,  76,  84,  87,
  22L,        "Developed",         62,  76,  84,  87,
  23L,        "Developed",         62,  76,  84,  87,
  24L,        "Developed",         62,  76,  84,  87,
  
  # Barren
  31L,        "Barren",            78,  86,  91,  93,
  
  # Forest merged: average of 41, 42, 43
  41L,        "Forest",            35,  60,  73,  80,
  42L,        "Forest",            35,  60,  73,  80,
  43L,        "Forest",            35,  60,  73,  80,
  
  # Shrub/scrub
  52L,        "Shrub/scrub",       55,  72,  81,  86,
  
  # Grassland
  71L,        "Grassland",         50,  69,  79,  84,
  
  # Agriculture merged: average of 81, 82
  81L,        "Agriculture",       58,  74,  82,  87,
  82L,        "Agriculture",       58,  74,  82,  87,
  
  # Wetlands merged: same original values for 90 and 95
  90L,        "Wetlands",          30,  58,  71,  78,
  95L,        "Wetlands",          30,  58,  71,  78
)
# Choose one assumed soil group for all pixels
assumed_hsg <- "C"

# Reclassification matrix: NLCD code -> CN
rcl_matrix <- cn_lookup %>%
  select(nlcd_code, CN = all_of(assumed_hsg)) %>%
  as.matrix()

# ----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ------------------------------------------------------------
# 1. Make group-level CN lookup using assumed HSG
# ------------------------------------------------------------

assumed_hsg <- "C"

cn_group_lookup <- cn_lookup %>%
  select(nlcd_group, CN = all_of(assumed_hsg)) %>%
  distinct(nlcd_group, CN) %>%
  mutate(
    group = case_when(
      nlcd_group == "Water"       ~ "water",
      nlcd_group == "Developed"   ~ "developed",
      nlcd_group == "Barren"      ~ "barren",
      nlcd_group == "Forest"      ~ "forest",
      nlcd_group == "Shrub/scrub" ~ "shrub_grass",
      nlcd_group == "Grassland"   ~ "shrub_grass",
      nlcd_group == "Agriculture" ~ "agriculture",
      nlcd_group == "Wetlands"    ~ "wetland"
    )
  ) %>%
  group_by(group) %>%
  summarise(CN = mean(CN), .groups = "drop")

# ------------------------------------------------------------
# 2. Convert lc_shares to long format and calculate HUC12 CN
# ------------------------------------------------------------
cn_huc12 <- lc_shares %>%
  as_tibble() %>%
  select(huc12, huc4, matches("_share_\\d{4}$")) %>%
  pivot_longer(
    cols = matches("_share_\\d{4}$"),
    names_to = c("group", "year"),
    names_pattern = "(.+)_share_(\\d{4})",
    values_to = "share"
  ) %>%
  mutate(
    year = as.integer(year),
    huc12 = str_pad(as.character(huc12), 12, pad = "0")
  ) %>%
  left_join(cn_group_lookup, by = "group") %>%
  group_by(huc12, huc4, year) %>%
  summarise(
    CN = sum(share * CN, na.rm = TRUE),
    share_check = sum(share, na.rm = TRUE),
    .groups = "drop"
  )

head(cn_huc12)

write_csv(cn_huc12, "data/cn/cn_huc12_sac.csv")


# --------------------------------------------------------------------------

# keep the basin-wide event definition, but compute within-basin spatial rainfall exposure only on those event days.
# 1. Load basin-wide event days
target_huc4s = "1802"
# run build_events.R
event_days <- readRDS("sac/data/eventdays_post2009_sac.rds") %>%
  as.data.table()

event_days[, date := as.Date(date)]

event_dates <- unique(event_days[, .(
  date, event_id, event_start, event_end, peak_date,
  extreme_event, is_peak_day
)])

# 2. Load HUC12 precipitation
precip_all <- readRDS("sac/data/precip_huc12s_sac.rds") %>%
  as.data.table()

precip_all[, date := as.Date(date)]
precip_all[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
precip_all[, huc4 := str_sub(huc12, 1, 4)]

precip_all <- precip_all[huc4 %in% target_huc4s]

# 3. Compute HUC12-specific rolling 3-day precip
precip_huc12 <- precip_all[
  order(huc12, date)
][
  , roll3_huc12 := zoo::rollsum(precip_mm, k = 3, fill = NA, align = "right"),
  by = huc12
]

# 4. Keep only basin-wide event dates
precip_event_huc12 <- event_dates[
  precip_huc12,
  on = "date",
  nomatch = 0
]

# 5. Assign NLCD/CN year
nlcd_years <- c(2009, 2012, 2015, 2018, 2021, 2024)

assign_lc_year <- function(d) {
  yr <- as.integer(format(d, "%Y"))
  vapply(
    yr,
    function(y) as.integer(max(nlcd_years[nlcd_years <= y])),
    integer(1L)
  )
}
precip_event_huc12[, nlcd_year := assign_lc_year(date)]
precip_event_huc12[, precip_3day_in := roll3_huc12 / 25.4]

# 6. Join HUC12 CN
cn_huc12_dt <- as.data.table(cn_huc12)
setnames(cn_huc12_dt, "year", "nlcd_year")

precip_event_huc12 <- merge(
  precip_event_huc12,
  cn_huc12_dt[, .(huc12, nlcd_year, CN)],
  by = c("huc12", "nlcd_year"),
  all.x = TRUE
)

# 7. Compute runoff
scs_runoff <- function(P, CN) {
  S  <- (1000 / CN) - 10
  Ia <- 0.2 * S
  
  fifelse(
    is.na(P) | is.na(CN),
    NA_real_,
    fifelse(
      P > Ia,
      (P - Ia)^2 / ((P - Ia) + S),
      0
    )
  )
}

precip_event_huc12[, runoff_in := scs_runoff(precip_3day_in, CN)]

# saveRDS(
#   precip_event_huc12,
#   "sac/data/runoff_huc12_sac.rds"
# )
