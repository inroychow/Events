setwd("L:\\Wetland Flood Mitigation\\Paper_NFIP")
# Load packages
pacman::p_load(sf, terra, exactextractr, data.table, dplyr, lubridate, stringr, prism, fixest, future, future.apply, progressr, tidyverse, fst, patchwork, scales)

#Load in base data
htdt = fread("data/htdt_allconus.xls")

nfip = fread("L:\\Wetland Flood Mitigation\\NFIP\\Claims\\FimaNfipClaimsV2.csv")
setDT(nfip)  # raw NFIP table

nfip_panel <- {
  nfip[, `:=`(
    date     = as.Date(dateOfLoss),
    tract_id = as.character(censusTract),
    loss_amt = rowSums(.SD, na.rm = TRUE)
  ), .SDcols = c("netBuildingPaymentAmount",
                 "netContentsPaymentAmount")][
                   , .(loss_amt = sum(loss_amt), n_claims = .N), by = .(tract_id, date)
                 ][
                   CJ(tract_id = unique(tract_id),
                      date     = seq(min(date), max(date), by = "day")),
                   on = .(tract_id, date)
                 ][
                   is.na(loss_amt), `:=`(loss_amt = 0, n_claims = 0)
                 ][
                   , any_claim := as.integer(n_claims > 0)
                 ][]
}



#### EXTRACT PRISM #####

# Choose a local folder 
prism_set_dl_dir("PRISM/prism_archive")  

# Date range to match nfip panel
min_d <- min(nfip_panel$date)
max_d <- max(nfip_panel$date)

# PRISM daily starts 1981-01-01
min_d <- max(min_d, as.Date("2009-01-01"))

# Download daily precipitation #4km resolution
get_prism_dailys(
  type     = "ppt",
  minDate  = min_d,
  maxDate  = max_d,
  keepZip  = FALSE,
  resolution = "4km"
)