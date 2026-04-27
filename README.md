
We investigate whether changes in upstream land use alter the relationship between heavy precipitation events and downstream flood damages. The core intuition is that the same rainfall event may generate different downstream flood outcomes depending on the composition of upstream land cover and how that land cover has changed over time. As a proof of concept, we focus initially on the Sacramento River Basin HUC4 watershed, using daily National Flood Insurance Program (NFIP) flood insurance claims from 2009 onward and daily precipitation from PRISM. We ultimately plan to scale this framework to all 250 HUC4 watersheds in the contiguous United States (CONUS).## 1. System Requirements

### Software Dependencies

We use **RStudio version 4.4.3** to run all code. Required R packages are loaded within the scripts. This code is run on the DEEP HPC server of UC Davis.

## 2. Installation Guide

### Instructions

1. Install **R**: https://cran.r-project.org/  
2. Install **RStudio**: https://posit.co/download/rstudio-desktop/ (typical install time is 5-10 minutes).
3. Clone or download this repository  
4. Open the project in RStudio
5. Install required packages, indicated in the individual script files
6. Set working directory to the project folder.

## Files and variables

### 3. Data
- **PRISM:** Daily precipitation data (2001–2024)
- **NLCD:** Land cover snapshots (2009, 2012, 2016, 2019, 2021, 2024)
- **NHDPlus:** Hydrologic network used to map upstream and downstream relationships between spatial units
- **NFIP:** Census tract-level flood insurance claims and policy data used as a proxy for flood damages.

## 4. Code and software
These files should be run in order. 

- **build_htdt_distances.R:** calculates km distance between upstream HUC12 centroids and downstream census tract centroids within HUC4 watersheds.
- **extract_prism.R:** extracts and saves daily precipitation data for the HUC4 watershed.
- **crop_prism.R:** processes PRISM daily precipitation rasters into a panel of mean daily precipitation at the HUC12 level for the selected HUC4 watershed.
-  **build_events.R:** builds extreme events (97th percentile) for the HUC4 watersheds, starting in the year 2009.
- **build_nfip.R:** cleans NFIP claims and policy data, aggregating to the census tract level.
- **build_lc .R:** builds land cover snapshots for the panel years for the HUC4 watershed, starting in the year 2009.
- **build_runoff:** uses USDA Soil Conservation Service Curve Numbers (CN) to estimate downstream runoff for each HUC12, based on land cover proportion.
- **build_weights:** builds weights based on NFIP policies in 2009, exponential distance decay, and HUC12 watershed area.
- **apply_weights:** applies weights from build_weights.R to land cover shares and runoff upstream of each census tract.
- **create_panel.R:** Generates main dataset for analysis, with land cover at the HUC12 watershed level, in event windows only.
- **build_huc4_panel.R:** Aggregates tract-level NFIP claims, policies, upstream runoff, land cover, and HUC12 precipitation into a HUC4-day event-study panel centered on extreme precipitation events, with precipitation weighted by baseline 2009 NFIP policy exposure.
- **analysis_sac.R:** Runs event study models of relationship between upstream land cover change and flood damages. 
  
