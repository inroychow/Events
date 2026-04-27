setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
library(sf)
library(terra)
library(data.table)
library(stringr)
library(dplyr)

# Need to run build_htdt_watershed.R first.

# target_huc4s <- c("1806", "1801", "1805", "1804", "1803")

target_huc4s = "1802"
# -----------------------------
# Load data
# -----------------------------
policies <- fread("data/NFIP/policies.csv")
claims   <- fread("data/NFIP/claims.csv")
htdt     <- readRDS("sac/data/htdt_1802.rds")

target_tracts = unique(htdt$tract)
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


