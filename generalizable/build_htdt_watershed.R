library(sf)
library(dplyr)
library(data.table)
library(stringr)

# -----------------------------
# Function: get downstream tracts within any HUC4 boundary
# -----------------------------

get_huc4_downstream_tracts <- function(target_huc4s,
                                       tracts_path = "data/tracts_2024.shp",
                                       wbd_path = "data/WBD_National_GDB.gdb",
                                       htdt) {
  
  target_huc4s <- str_pad(as.character(target_huc4s), width = 4, pad = "0")
  
  # Load spatial boundaries
  tracts_sf <- st_read(tracts_path, quiet = TRUE)
  watersheds <- st_read(wbd_path, layer = "WBDHU12", quiet = TRUE)
  
  # Make sure HUC12 field exists and is padded
  watersheds <- watersheds %>%
    mutate(
      huc12 = str_pad(as.character(huc12), width = 12, pad = "0"),
      huc4 = substr(huc12, 1, 4)
    )
  
  # Dissolve selected HUC12s into HUC4 boundaries
  huc4_target <- watersheds %>%
    filter(huc4 %in% target_huc4s) %>%
    group_by(huc4) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
  
  # Transform tracts to HUC CRS
  tracts_sf <- st_transform(tracts_sf, st_crs(huc4_target))
  
  # Use representative points instead of centroids to avoid points outside polygons
  tract_points <- st_point_on_surface(tracts_sf)
  
  tracts_in_huc4 <- st_join(
    tract_points,
    huc4_target["huc4"],
    join = st_within
  ) %>%
    filter(!is.na(huc4)) %>%
    st_drop_geometry() %>%
    transmute(
      tract = str_pad(as.character(GEOID), width = 11, pad = "0"),
      huc4 = str_pad(as.character(huc4), width = 4, pad = "0")
    ) %>%
    distinct()
  
  setDT(tracts_in_huc4)
  
  # -----------------------------
  # Build HTDT long crosswalk
  # -----------------------------
  
  htdt <- as.data.table(htdt)
  
  htdt_long <- htdt[, .(
    tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ",\\s*"))
  ), by = huc12]
  
  htdt_long[, `:=`(
    tract = str_pad(as.character(tract), width = 11, pad = "0"),
    huc12 = str_pad(as.character(huc12), width = 12, pad = "0")
  )]
  
  htdt_long[, huc4 := substr(huc12, 1, 4)]
  
  # Keep only requested HUC4s
  htdt_target <- htdt_long[huc4 %in% target_huc4s]
  
  # Keep only downstream tracts that are spatially inside the same HUC4
  tract_huc4 <- merge(
    unique(htdt_target[, .(tract, huc12, huc4)]),
    tracts_in_huc4,
    by = c("tract", "huc4")
  )
  
  target_tracts <- unique(tract_huc4$tract)
  
  return(list(
    tract_huc4 = tract_huc4,
    target_tracts = target_tracts,
    huc4_boundaries = huc4_target
  ))
}


result <- get_huc4_downstream_tracts(
  target_huc4s = "1802",
  htdt = htdt
)

saveRDS(result, "sac/data/htdt_1802.rds")