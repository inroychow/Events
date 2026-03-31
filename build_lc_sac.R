
setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")

library(sf)          
library(terra)       
library(data.table)  
library(stringr)     
library(dplyr)       

watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)
huc12_sac <- watersheds[str_starts(watersheds$huc12, "1802"), ]

geom <- st_geometry(huc12_sac)
huc12_sac <- st_drop_geometry(huc12_sac) %>%
  mutate(across(where(is.list), ~ vapply(., \(x) paste(unlist(x), collapse = "|"), character(1))))
st_geometry(huc12_sac) <- geom

huc12_sac <- st_as_sf(huc12_sac) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4269)

cat(sprintf("Sacramento HUC12s: %d polygons\n", nrow(huc12_sac)))

nlcd_groups <- list(
  water = 11,
  developed = c(21, 22, 23, 24),
  barren = 31,
  forest = c(41, 42, 43),
  shrub_grass = c(51, 52, 71, 72, 73, 74),
  agriculture = c(81, 82),
  wetland = c(90, 95)
)

compute_lc_shares <- function(nlcd_rast, huc12_rast, groups, yr) {
  pixel_area_sqkm <- prod(res(nlcd_rast)) / 1e6
  
  total_cells <- zonal(!is.na(nlcd_rast), huc12_rast, fun = "sum", na.rm = TRUE) |>
    as.data.table()
  setnames(total_cells, c("huc12", "total_cells"))
  
  group_dt <- rbindlist(lapply(names(groups), \(grp) {
    out <- zonal(as.numeric(nlcd_rast %in% groups[[grp]]), huc12_rast, fun = "sum", na.rm = TRUE) |>
      as.data.table()
    setnames(out, c("huc12", "n_cells"))
    out[, group := grp]
    out
  }))
  
  group_dt <- merge(group_dt, total_cells, by = "huc12")
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
  
  merge(share_wide, area_wide, by = "huc12")
}

# years 
years <- c(2009, 2012, 2015, 2018, 2021, 2024)

# use first year raster as template
template_rast <- rast(sprintf("NLCD rasters/nlcd%d.tif", years[1]))
huc12_vect <- project(vect(huc12_sac), crs(template_rast))

template_down <- aggregate(
  mask(crop(template_rast, huc12_vect), huc12_vect),
  fact = 3,
  fun = "modal"
)

huc12_rast <- rasterize(huc12_vect, template_down, field = "huc12")

# compute land cover shares for each year
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

# merge all years together
lc_shares <- Reduce(function(x, y) merge(x, y, by = "huc12", all = TRUE), lc_list)


saveRDS(lc_shares, "data/lc_shares_sac.rds")

# ----------------------------------

# LC shares by tract sac---------------------

# -------------------------------------

lc_shares = readRDS("data/lc_shares_sac.rds")


htdt <- fread("data/htdt_allconus.xls")
htdt_long <- htdt[, .(tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ", "))), by = huc12]
htdt_long[, tract := as.character(tract)]
htdt_sac <- htdt_long %>% filter(huc12 %in% huc12_sac$huc12)
htdt_sac[, huc12 := str_pad(as.character(huc12), width = 12, pad = "0")]
htdt_sac <- htdt_sac %>% rename(upstream_huc12 = huc12)

# ---


lc_shares$huc12 <- str_pad(as.character(lc_shares$huc12), 12, pad = "0")
htdt_sac$upstream_huc12 <- str_pad(as.character(htdt_sac$upstream_huc12), 12, pad = "0")
htdt_sac$tract <- as.character(htdt_sac$tract)

huc_area <- huc12_sac %>%
  st_transform(5070) %>%
  mutate(area_sqkm_huc = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(huc12, area_sqkm_huc)

lc_long <- lc_shares %>%
  pivot_longer(
    cols = -huc12,
    names_to = c("group", "metric", "year"),
    names_pattern = "(.+)_(share|sqkm)_(\\d{4})",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = c(group, metric), values_from = value)

tract_huc_lc <- htdt_sac %>%
  left_join(lc_long, by = c("upstream_huc12" = "huc12")) %>%
  left_join(huc_area, by = c("upstream_huc12" = "huc12"))

value_cols <- tract_huc_lc %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff("area_sqkm_huc")

tract_year_lc <- tract_huc_lc %>%
  group_by(tract, year) %>%
  summarise(
    across(all_of(value_cols), ~ weighted.mean(.x, w = area_sqkm_huc, na.rm = TRUE)),
    .groups = "drop"
  )

# saveRDS(tract_year_lc, "data/tract_year_lc.rds")


# --------------------------------------------------
# Simplified: developed vs. undeveloped only
# NLCD codes 21-24 = developed; everything else = undeveloped
# --------------------------------------------------
lc_shares <- readRDS("data/lc_shares_sac.rds")
lc_shares$huc12 <- str_pad(as.character(lc_shares$huc12), 12, pad = "0")

# Derive undeveloped_share as complement of developed_share, per year
for (yr in c(2009, 2012, 2015, 2018, 2021, 2024)) {
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
# Pivot to long
lc_long <- lc_shares %>%
  select(huc12, matches("^(developed|undeveloped)_(share|sqkm)_\\d{4}$")) %>%
  pivot_longer(
    cols = -huc12,
    names_to  = c("group", "metric", "year"),
    names_pattern = "(.+)_(share|sqkm)_(\\d{4})",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = c(group, metric), values_from = value)

# Area-weighted average upstream land cover per tract-year
tract_year_lc <- htdt_sac %>%
  left_join(lc_long,  by = c("upstream_huc12" = "huc12")) %>%
  left_join(huc_area, by = c("upstream_huc12" = "huc12")) %>%
  group_by(tract, year) %>%
  summarise(
    across(c(developed_share, undeveloped_share, developed_sqkm, undeveloped_sqkm),
           ~ weighted.mean(.x, w = area_sqkm_huc, na.rm = TRUE)),
    .groups = "drop"
  )

saveRDS(tract_year_lc, "sac/data/tract_year_lc2.rds")

