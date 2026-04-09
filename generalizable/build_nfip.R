setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
library(sf)
library(terra)
library(data.table)
library(stringr)
library(dplyr)

# target_huc4s <- c("1806", "1801", "1805", "1804", "1803")

target_huc4s = "1802"
# -----------------------------
# Load data
# -----------------------------
policies <- fread("data/NFIP/policies.csv")
claims   <- fread("data/NFIP/claims.csv")
htdt     <- fread("data/htdt_allconus.xls")

# Load spatial boundaries
tracts_sf <- st_read("data/tracts_2024.shp")   
watersheds <- st_read("data/WBD_National_GDB.gdb", layer = "WBDHU12", quiet = TRUE)

# Derive HUC4 boundaries by dissolving HUC12s on their first 4 characters
watersheds[, "huc4"] <- substr(watersheds$huc12, 1, 4) 

huc4_target <- watersheds %>%
  filter(huc4 %in% target_huc4s) %>%
  group_by(huc4) %>%
  summarise(geometry = st_union(shape), .groups = "drop")
# -----------------------------
# Build spatially-bounded tract-HUC4 crosswalk
# -----------------------------

# Spatially join tracts to HUC4s — only keep tracts whose centroid/area
#    falls within a target HUC4
tracts_sf <- st_transform(tracts_sf, st_crs(huc4_target))

# Use centroid join to avoid edge-overlap double-counting
tract_centroids <- st_centroid(tracts_sf)
tracts_in_huc4 <- st_join(tract_centroids, huc4_target["huc4"], join = st_within)
tracts_in_huc4 <- tracts_in_huc4 %>%
  filter(!is.na(huc4)) %>%
  st_drop_geometry() %>%
  select(tract = GEOID, huc4) %>%   # adjust GEOID column name as needed
  distinct()

setDT(tracts_in_huc4)
tracts_in_huc4[, `:=`(
  tract = str_pad(as.character(tract), width = 11, pad = "0"),
  huc4  = str_pad(as.character(huc4),  width = 4,  pad = "0")
)]

# -----------------------------
# Build HTDT crosswalk — but NOW intersect with spatially-bounded tracts
# -----------------------------
htdt_long <- htdt[, .(
  tract = unlist(strsplit(gsub("\\[|\\]|'", "", downstream_tracts), ", "))
), by = huc12]

htdt_long[, `:=`(
  tract = str_pad(as.character(tract), width = 11, pad = "0"),
  huc12 = str_pad(as.character(huc12), width = 12, pad = "0")
)]
htdt_long[, huc4 := substr(huc12, 1, 4)]

# Keep only target HUC4s in HTDT
htdt_target <- htdt_long[huc4 %in% target_huc4s]

# KEY FIX: inner join with spatially-bounded tracts
# This drops any downstream tract that lies outside the HUC4 boundary
tract_huc4 <- merge(
  unique(htdt_target[, .(tract, huc4)]),
  tracts_in_huc4,
  by = c("tract", "huc4")
)

target_tracts <- unique(tract_huc4$tract)

# -----------------------------
# Claims
# -----------------------------
claims[, tract := str_pad(as.character(censusTract), width = 11, pad = "0")]
claims[, dateOfLoss := as.Date(dateOfLoss)]

claims_target <- claims[tract %in% target_tracts]
claims_target <- merge(
  claims_target,
  tract_huc4,
  by = "tract",
  all.x = FALSE, all.y = FALSE,
  allow.cartesian = TRUE
)

claims_target <- claims_target[, .(
  tract_claims     = .N,
  tract_total_paid = sum(amountPaidOnBuildingClaim + amountPaidOnContentsClaim, na.rm = TRUE)
), by = .(tract, huc4, dateOfLoss)]

# -----------------------------
# Policies
# -----------------------------
policies <- policies[, .(
  censusTract,
  policyEffectiveDate,
  policyTerminationDate,
  total_coverage = totalBuildingInsuranceCoverage + totalContentsInsuranceCoverage
)]
policies[, censusTract := str_pad(as.character(censusTract), width = 11, pad = "0")]
setnames(policies, "censusTract", "tract")
policies[, `:=`(
  policyEffectiveDate    = as.Date(policyEffectiveDate),
  policyTerminationDate  = as.Date(policyTerminationDate)
)]

policies_target <- policies[tract %in% target_tracts]
policies_target <- merge(
  policies_target,
  tract_huc4,
  by = "tract",
  all.x = FALSE, all.y = FALSE,
  allow.cartesian = TRUE
)

policies_target[, n_policies_tract_huc4 := .N, by = .(tract, huc4)]

# -----------------------------
# Save
# -----------------------------
# saveRDS(claims_target,  "sac/data/claims_sac.rds")
# saveRDS(policies_target, "sac/data/policies_sac.rds")


