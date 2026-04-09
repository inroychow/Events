setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")

library(sf)
library(terra)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)

# -----------------------------
# User inputs
# -----------------------------
# target_huc4s <- c("1806", "1801", "1805", "1804", "1803")
target_huc4s <- "1802"
years <- c(2009, 2012, 2015, 2018, 2021, 2024)

nlcd_groups <- list(
  water = 11,
  developed = c(21, 22, 23, 24),
  barren = 31,
  forest = c(41, 42, 43),
  shrub_grass = c(51, 52, 71, 72, 73, 74),
  agriculture = c(81, 82),
  wetland = c(90, 95)
)

output_stub <- paste(target_huc4s, collapse = "_")

# -----------------------------
# Helper function
# -----------------------------
compute_lc_shares <- function(nlcd_rast, huc12_rast, groups, yr) {
  pixel_area_sqkm <- prod(res(nlcd_rast)) / 1e6
  
  total_cells <- zonal(!is.na(nlcd_rast), huc12_rast, fun = "sum", na.rm = TRUE) |>
    as.data.table()
  setnames(total_cells, c("huc12", "total_cells"))
  
  group_dt <- rbindlist(lapply(names(groups), function(grp) {
    out <- zonal(as.numeric(nlcd_rast %in% groups[[grp]]), huc12_rast, fun = "sum", na.rm = TRUE) |>
      as.data.table()
    setnames(out, c("huc12", "n_cells"))
    out[, group := grp]
    out
  }))
  
  group_dt <- merge(group_dt, total_cells, by = "huc12", all.x = TRUE)
  
  group_dt[, `:=`(
    share = n_cells / total_cells,
    area_sqkm = n_cells * pixel_area_sqkm
  )]
  
  share_wide <- dcast(group_dt, huc12 ~ group, value.var = "share")
  area_wide  <- dcast(group_dt, huc12 ~ group, value.var = "area_sqkm")
  
  setnames(
    share_wide,
    setdiff(names(share_wide), "huc12"),
    paste0(setdiff(names(share_wide), "huc12"), "_share_", yr)
  )
  
  setnames(
    area_wide,
    setdiff(names(area_wide), "huc12"),
    paste0(setdiff(names(area_wide), "huc12"), "_sqkm_", yr)
  )
  
  merge(share_wide, area_wide, by = "huc12", all = TRUE)
}

# -----------------------------
# Read and prep HUC12 polygons
# -----------------------------
watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

huc12_subset <- watersheds[str_sub(watersheds$huc12, 1, 4) %in% target_huc4s, ]

geom <- st_geometry(huc12_subset)
huc12_subset <- st_drop_geometry(huc12_subset) %>%
  mutate(across(where(is.list), ~ vapply(., function(x) paste(unlist(x), collapse = "|"), character(1))))
st_geometry(huc12_subset) <- geom

huc12_subset <- st_as_sf(huc12_subset) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4269) %>%
  mutate(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    huc4  = str_sub(huc12, 1, 4)
  )

cat(sprintf("Selected HUC12s: %d polygons across %d HUC4s\n",
            nrow(huc12_subset), length(unique(huc12_subset$huc4))))

# -----------------------------
# Read tract boundaries and keep only tracts in target HUC4s
# -----------------------------
tracts_sf <- st_read("data/tracts_2024.shp", quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(st_crs(huc12_subset))

# CHANGE GEOID IF YOUR TRACT ID FIELD HAS A DIFFERENT NAME
tracts_sf <- tracts_sf %>%
  mutate(tract = as.character(GEOID))

# dissolve HUC12s to HUC4 polygons
huc4_target <- huc12_subset %>%
  select(huc4) %>%
  group_by(huc4) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_make_valid()

tracts_sf <- st_make_valid(tracts_sf)

# assign each tract to a HUC4 using centroid
# this avoids problems if a tract touches multiple HUC4s
tract_centroids <- st_centroid(tracts_sf)

tract_huc4 <- st_join(
  tract_centroids %>% select(tract),
  huc4_target %>% select(huc4),
  join = st_within,
  left = FALSE
) %>%
  st_drop_geometry() %>%
  distinct(tract, huc4)

# optional: keep only tract polygons that belong to target HUC4s
tracts_sf_subset <- tracts_sf %>%
  inner_join(tract_huc4, by = "tract")

cat(sprintf("Selected tracts in target HUC4s: %d\n", nrow(tracts_sf_subset)))

# -----------------------------
# Build template raster
# -----------------------------
template_rast <- rast(sprintf("NLCD rasters/nlcd%d.tif", years[1]))
huc12_vect <- project(vect(huc12_subset), crs(template_rast))

template_down <- aggregate(
  mask(crop(template_rast, huc12_vect), huc12_vect),
  fact = 3,
  fun = "modal"
)

huc12_rast <- rasterize(huc12_vect, template_down, field = "huc12")

# -----------------------------
# Compute land cover by HUC12 for each year
# -----------------------------
lc_list <- lapply(years, function(yr) {
  cat("Processing year:", yr, "\n")
  
  nlcd <- rast(sprintf("NLCD rasters/nlcd%d.tif", yr))
  
  nlcd_down <- aggregate(
    mask(crop(nlcd, huc12_vect), huc12_vect),
    fact = 3,
    fun = "modal"
  )
  
  compute_lc_shares(nlcd_down, huc12_rast, nlcd_groups, yr)
})

lc_shares <- Reduce(function(x, y) merge(x, y, by = "huc12", all = TRUE), lc_list) %>%
  as.data.table()

lc_shares[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
lc_shares[, huc4 := str_sub(huc12, 1, 4)]

# -----------------------------
# Add HUC12 polygon area
# -----------------------------
huc_area <- huc12_subset %>%
  st_transform(5070) %>%
  mutate(area_sqkm_huc12 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  transmute(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    area_sqkm_huc12
  )

lc_shares <- merge(lc_shares, huc_area, by = "huc12", all.x = TRUE)

# -----------------------------
# Add undeveloped variables
# -----------------------------
for (yr in years) {
  lc_shares[[paste0("undeveloped_share_", yr)]] <-
    1 - lc_shares[[paste0("developed_share_", yr)]]
  
  lc_shares[[paste0("undeveloped_sqkm_", yr)]] <-
    lc_shares[[paste0("water_sqkm_", yr)]] +
    lc_shares[[paste0("barren_sqkm_", yr)]] +
    lc_shares[[paste0("forest_sqkm_", yr)]] +
    lc_shares[[paste0("shrub_grass_sqkm_", yr)]] +
    lc_shares[[paste0("agriculture_sqkm_", yr)]] +
    lc_shares[[paste0("wetland_sqkm_", yr)]]
}

head(lc_shares)
# saveRDS(lc_shares, "sac/data/lc_huc12_sac.rds")
# -----------------------------
# Build upstream HUC12 -> tract table
# and filter so tract is in same HUC4 as upstream HUC12
# -----------------------------
htdt <- fread("data/htdt_allconus.xls")

htdt_long <- htdt[, .(
  tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ", "))
), by = huc12]

htdt_long[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
htdt_long[, tract := as.character(trimws(tract))]

htdt_subset <- htdt_long[huc12 %in% huc12_subset$huc12]
setnames(htdt_subset, "huc12", "upstream_huc12")

# add huc4 from upstream huc12
htdt_subset[, huc4 := str_sub(upstream_huc12, 1, 4)]

# keep only rows where tract belongs to that same huc4
htdt_subset <- htdt_subset %>%
  inner_join(tract_huc4, by = c("tract", "huc4"))

cat(sprintf("Rows in filtered HUC12-tract table: %d\n", nrow(htdt_subset)))
cat(sprintf("Unique filtered tracts: %d\n", dplyr::n_distinct(htdt_subset$tract)))

# area table for weighting
huc_area_weight <- huc12_subset %>%
  st_transform(5070) %>%
  mutate(area_sqkm_huc = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  transmute(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    area_sqkm_huc
  )

# -----------------------------
# Pivot HUC12 LC to long by year
# -----------------------------
lc_long <- lc_shares %>%
  pivot_longer(
    cols = matches("_(share|sqkm)_\\d{4}$"),
    names_to = c("group", "metric", "year"),
    names_pattern = "(.+)_(share|sqkm)_(\\d{4})",
    values_to = "value"
  ) %>%
  pivot_wider(
    id_cols = c(huc12, huc4, year),
    names_from = c(group, metric),
    values_from = value
  ) %>%
  mutate(year = as.integer(year))

# -----------------------------
# Join upstream HUC12 LC to tracts
# -----------------------------
tract_huc_lc <- htdt_subset %>%
  left_join(lc_long, by = c("upstream_huc12" = "huc12", "huc4" = "huc4")) %>%
  left_join(huc_area_weight, by = c("upstream_huc12" = "huc12"))

value_cols <- tract_huc_lc %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(c("area_sqkm_huc", "year"))

# -----------------------------
# Tract-year summaries
# -----------------------------
tract_year_lc_areaweighted <- tract_huc_lc %>%
  group_by(tract, year) %>%
  summarise(
    across(all_of(value_cols), ~ weighted.mean(.x, w = area_sqkm_huc, na.rm = TRUE)),
    .groups = "drop"
  )

tract_year_lc_unweighted <- tract_huc_lc %>%
  group_by(tract, year) %>%
  summarise(
    across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
# 
# saveRDS(
#   tract_year_lc_areaweighted,
#   "sac/data/tract_year_lc_sac_aw.rds"
# )
# 
# saveRDS(
#   tract_year_lc_unweighted,
#   "sac/data/tract_year_lc_sac.rds"
# )


