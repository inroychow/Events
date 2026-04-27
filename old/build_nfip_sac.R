
setwd("L:\\Wetland Flood Mitigation\\Paper_claims")

library(sf)          
library(terra)       
library(data.table)  
library(stringr)     
library(dplyr)   

# Load data 

policies <- fread("data/policies.csv")
claims = fread("data/claims/claims.csv")
htdt = fread("data/htdt_allconus.xls")
# claims 
claims[, tract := str_pad(as.character(censusTract), width = 11, pad = "0")]
claims_sac <- claims %>% filter(tract %in% htdt_sac$tract)
claims_sac <- claims_sac[, .(
  tract_claims = .N,
  tract_total_paid = sum(amountPaidOnBuildingClaim + amountPaidOnContentsClaim, na.rm = TRUE)
), by = .(tract, dateOfLoss)]

#policies
policies <- policies[, .(
  censusTract, policyEffectiveDate, policyTerminationDate,
  total_coverage = totalBuildingInsuranceCoverage + totalContentsInsuranceCoverage
)]
policies[, censusTract := str_pad(as.character(censusTract), width = 11, pad = "0")]
policies[, n_policies_tract := .N, by = censusTract]
policies = policies %>% 
  rename(tract = censusTract)
policies_sac <- policies %>% filter(tract %in% htdt_sac$tract)

saveRDS(claims_sac, "sac/data/claims_sac.rds")
saveRDS(policies_sac, "sac/data/policies_sac.rds")

