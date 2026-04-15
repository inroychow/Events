setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")

library(sf)
library(data.table)
library(dplyr)
library(stringr)

# -----------------------------
# User inputs
# -----------------------------
out_dir <- "data/htdt_distance_by_huc4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read and prep HUC12 polygons
# -----------------------------
watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

geom <- st_geometry(watersheds)
watersheds <- st_drop_geometry(watersheds) %>%
  mutate(across(where(is.list), ~ vapply(., function(x) paste(unlist(x), collapse = "|"), character(1))))
st_geometry(watersheds) <- geom

huc12_sf <- st_as_sf(watersheds) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4269) %>%
  mutate(
    huc12 = str_pad(as.character(huc12), 12, pad = "0"),
    huc4  = str_sub(huc12, 1, 4)
  ) %>%
  select(huc12, huc4) %>% 
  st_make_valid

cat(sprintf("Loaded %d HUC12 polygons across %d HUC4s\n",
            nrow(huc12_sf), dplyr::n_distinct(huc12_sf$huc4)))

# -----------------------------
# Read and prep tract polygons
# -----------------------------
tracts_sf <- st_read("data/tracts_2024.shp", quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(st_crs(huc12_sf)) %>%
  mutate(tract = as.character(GEOID)) %>%
  st_make_valid() %>%
  select(tract)

cat(sprintf("Loaded %d tract polygons\n", nrow(tracts_sf)))

# -----------------------------
# Assign each tract to a HUC4 using tract centroid
# -----------------------------
huc4_sf <- huc12_sf %>%
  select(huc4) %>%
  group_by(huc4) %>%
  summarise(do_union = TRUE, .groups = "drop") %>%
  st_make_valid()

tract_centroids_for_join <- st_centroid(tracts_sf)

tract_huc4 <- st_join(
  tract_centroids_for_join %>% select(tract),
  huc4_sf %>% select(huc4),
  join = st_within,
  left = FALSE
) %>%
  st_drop_geometry() %>%
  distinct(tract, huc4)

cat(sprintf("Assigned %d tracts to HUC4s\n", nrow(tract_huc4)))

# -----------------------------
# Read and expand htdt
# -----------------------------
htdt <- fread("data/htdt_allconus.xls")

htdt_long <- htdt[, .(
  tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ", "))
), by = huc12]

htdt_long[, huc12 := str_pad(as.character(huc12), 12, pad = "0")]
htdt_long[, tract := as.character(trimws(tract))]
setnames(htdt_long, "huc12", "upstream_huc12")

# add huc4 from upstream huc12
htdt_long[, huc4 := str_sub(upstream_huc12, 1, 4)]

# keep only rows where tract belongs to that same huc4
htdt_long <- htdt_long %>%
  inner_join(tract_huc4, by = c("tract", "huc4")) %>%
  distinct(upstream_huc12, tract, huc4)

cat(sprintf("Filtered to %d valid upstream_huc12-tract pairs\n", nrow(htdt_long)))

# -----------------------------
# Build centroid layers in projected CRS for distance
# -----------------------------
huc12_centroids <- huc12_sf %>%
  st_transform(5070) %>%
  st_centroid() %>%
  select(upstream_huc12 = huc12, huc4)

tract_centroids <- tracts_sf %>%
  inner_join(tract_huc4, by = "tract") %>%
  st_transform(5070) %>%
  st_centroid() %>%
  select(tract, huc4)

# -----------------------------
# Loop over HUC4s and save one file per HUC4
# -----------------------------
all_huc4s <- sort(unique(htdt_long$huc4))

for (this_huc4 in all_huc4s) {
  cat("\n-----------------------------\n")
  cat("Processing HUC4:", this_huc4, "\n")
  
  pairs_h4 <- htdt_long %>%
    filter(huc4 == this_huc4)
  
  if (nrow(pairs_h4) == 0) {
    cat("No pairs found, skipping\n")
    next
  }
  
  huc12_h4 <- huc12_centroids %>%
    filter(huc4 == this_huc4) %>%
    rename(huc12_geom = geometry)
  
  tracts_h4 <- tract_centroids %>%
    filter(huc4 == this_huc4) %>%
    rename(tract_geom = geometry)
  
  # attach geometries as list-columns
  pairs_geom <- pairs_h4 %>%
    left_join(
      huc12_h4 %>%
        st_drop_geometry() %>%
        mutate(huc12_geom = st_geometry(huc12_h4)),
      by = c("upstream_huc12", "huc4")
    ) %>%
    left_join(
      tracts_h4 %>%
        st_drop_geometry() %>%
        mutate(tract_geom = st_geometry(tracts_h4)),
      by = c("tract", "huc4")
    )
  
  # keep only rows with both geometries present
  pairs_geom <- pairs_geom %>%
    filter(!sapply(huc12_geom, is.null), !sapply(tract_geom, is.null))
  
  if (nrow(pairs_geom) == 0) {
    cat("No matched geometries, skipping\n")
    next
  }
  
  # convert geometry list-columns to sfc vectors
  huc12_sfc <- st_sfc(pairs_geom$huc12_geom, crs = 5070)
  tract_sfc <- st_sfc(pairs_geom$tract_geom, crs = 5070)
  
  # compute centroid-to-centroid distance in km
  pairs_geom$dist_km <- as.numeric(
    st_distance(huc12_sfc, tract_sfc, by_element = TRUE)
  ) / 1000
  
  dist_out <- pairs_geom %>%
    select(upstream_huc12, tract, huc4, dist_km) %>%
    arrange(upstream_huc12, tract)
  
  out_file <- file.path(out_dir, paste0("htdt_dist_", this_huc4, ".csv"))
  fwrite(dist_out, out_file)
  
  cat(sprintf("Saved %d rows to %s\n", nrow(dist_out), out_file))
}

cat("\nDone.\n")


